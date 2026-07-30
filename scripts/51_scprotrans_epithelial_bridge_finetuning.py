#!/usr/bin/env python3

"""Exploratory epithelial-bridge fine-tuning with frozen external test arms."""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scvi
import torch
from scipy import sparse
from sklearn.linear_model import Ridge


ROOT = Path(__file__).resolve().parents[1]
TRAINING_SCRIPT = ROOT / "scripts" / "39_scprotrans_crossdomain_training.py"


def load_module():
    spec = importlib.util.spec_from_file_location("scprotrans_crossdomain", TRAINING_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {TRAINING_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--application-dir",
        default=str(
            ROOT / "outputs/scprotrans_hgsoc_v4/scprotrans_benchmark/pilot/"
            "h500_sv10_pt10_seed260716/application"
        ),
    )
    parser.add_argument(
        "--prepared-dir",
        default=str(
            ROOT / "outputs/scprotrans_hgsoc_v4/scprotrans_benchmark/prepared/"
            "g128_r250_h2600_s260716"
        ),
    )
    parser.add_argument(
        "--bridge-dir",
        default=str(
            ROOT / "outputs/scprotrans_hgsoc_v4/scprotrans_benchmark/gse253719_bridge"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=str(
            ROOT / "outputs/scprotrans_hgsoc_v4/scprotrans_benchmark/"
            "epithelial_bridge_finetuning"
        ),
    )
    parser.add_argument("--seeds", nargs="+", type=int, default=[260718, 260719, 260720, 260721, 260722])
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--learning-rate", type=float, default=1e-4)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--split-seed", type=int, default=260718)
    return parser.parse_args()


def dense(values) -> np.ndarray:
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values, dtype=np.float32)


def load_latent(path: Path) -> np.ndarray:
    model = scvi.model.SCVI.load(str(path), accelerator="cpu", device="auto")
    return np.asarray(model.get_latent_representation(batch_size=512), dtype=np.float32)


def split_bridge_donors(donors: np.ndarray, seed: int) -> dict[str, np.ndarray]:
    unique = np.unique(donors)
    if unique.size < 9:
        raise RuntimeError("Bridge fine-tuning requires at least nine donors")
    rng = np.random.default_rng(seed)
    shuffled = unique.copy()
    rng.shuffle(shuffled)
    test_donors = set(shuffled[:3])
    validation_donors = set(shuffled[3:6])
    train_donors = set(shuffled[6:])
    return {
        "train": np.isin(donors, sorted(train_donors)),
        "validation": np.isin(donors, sorted(validation_donors)),
        "test": np.isin(donors, sorted(test_donors)),
        "train_donors": np.asarray(sorted(train_donors)),
        "validation_donors": np.asarray(sorted(validation_donors)),
        "test_donors": np.asarray(sorted(test_donors)),
    }


def fine_tune_hla(
    base_model,
    train_expression,
    train_truth,
    train_latent,
    validation_expression,
    validation_truth,
    validation_latent,
    gene_embedding,
    protein_embedding,
    target_index: int,
    *,
    seed: int,
    epochs: int,
    patience: int,
    learning_rate: float,
    batch_size: int,
):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = copy.deepcopy(base_model).to("cpu")
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    criterion = torch.nn.MSELoss()
    dataset = torch.utils.data.TensorDataset(
        torch.from_numpy(dense(train_expression)),
        torch.from_numpy(dense(train_truth).reshape(-1)),
        torch.from_numpy(dense(train_latent)),
    )
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=min(batch_size, len(dataset)),
        shuffle=True,
        generator=torch.Generator().manual_seed(seed),
    )
    gene_tensor = torch.from_numpy(dense(gene_embedding))
    protein_tensor = torch.from_numpy(dense(protein_embedding))
    val_expression = torch.from_numpy(dense(validation_expression))
    val_truth = torch.from_numpy(dense(validation_truth).reshape(-1))
    val_latent = torch.from_numpy(dense(validation_latent))
    best_state = None
    best_loss = float("inf")
    stale = 0
    history = []
    for epoch in range(epochs):
        model.train()
        total_loss = 0.0
        total_cells = 0
        for expression, truth, latent in loader:
            optimizer.zero_grad(set_to_none=True)
            prediction, _ = model(protein_tensor, gene_tensor, latent, expression)
            loss = criterion(prediction[:, target_index], truth)
            loss.backward()
            optimizer.step()
            total_loss += float(loss.detach()) * expression.shape[0]
            total_cells += expression.shape[0]
        model.eval()
        with torch.no_grad():
            val_prediction, _ = model(
                protein_tensor, gene_tensor, val_latent, val_expression
            )
            validation_loss = float(
                criterion(val_prediction[:, target_index], val_truth).detach()
            )
        history.append(
            {
                "epoch": epoch + 1,
                "train_loss": total_loss / total_cells,
                "validation_loss": validation_loss,
            }
        )
        if validation_loss < best_loss - 1e-8:
            best_loss = validation_loss
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
            stale = 0
        else:
            stale += 1
            if stale >= patience:
                break
    if best_state is None:
        raise RuntimeError("Fine-tuning did not produce a finite model")
    model.load_state_dict(best_state, strict=True)
    model.eval()
    return model, history


def unit_metrics(module, truth, predictions, unit: str, seed: int) -> list[dict[str, object]]:
    rows = []
    for method, values in predictions.items():
        rows.append(
            {
                "seed": seed,
                "unit": unit,
                "method": method,
                "n_cells": len(truth),
                **module.compute_prediction_metrics(truth, values),
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    module = load_module()
    application_dir = Path(args.application_dir)
    prepared_dir = Path(args.prepared_dir)
    bridge_dir = Path(args.bridge_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_model, metadata = module.load_model_bundle(application_dir / "scprotrans_model_bundle.pt")
    gene_order = list(metadata["gene_order"])
    target_index = list(metadata["protein_order"]).index("HLA_DR_COMPLEX")
    gene_embedding = np.load(prepared_dir / "gene_embedding.npy")
    protein_embedding = np.load(prepared_dir / "protein_embedding.npy")

    reference = module.load_prepared_block(prepared_dir, "reference")
    gse128639 = module.load_prepared_block(prepared_dir, "gse128639")
    gse254985 = module.load_prepared_block(prepared_dir, "gse254985")
    hgsoc = module.load_prepared_block(prepared_dir, "hgsoc")
    bridge_archive = np.load(bridge_dir / "gse253719_bridge_arrays.npz")
    bridge = {key: np.asarray(bridge_archive[key]) for key in bridge_archive.files}
    split = split_bridge_donors(bridge["donors"].astype(str), args.split_seed)

    reference_latent = load_latent(application_dir / "scvi_reference")
    query_latent = {
        "gse128639": load_latent(application_dir / "scvi_gse128639"),
        "gse254985": load_latent(application_dir / "scvi_gse254985"),
        "hgsoc": load_latent(application_dir / "scvi_hgsoc"),
    }
    reference_target = dense(reference["protein"])[:, target_index]
    train_expression = np.vstack(
        [dense(reference["expression"]), dense(bridge["expression"])[split["train"]]]
    )
    train_truth = np.concatenate(
        [reference_target, dense(bridge["truth"])[split["train"]]]
    )
    train_latent = np.vstack(
        [reference_latent, dense(bridge["latent"])[split["train"]]]
    )
    validation_expression = dense(bridge["expression"])[split["validation"]]
    validation_truth = dense(bridge["truth"])[split["validation"]]
    validation_latent = dense(bridge["latent"])[split["validation"]]

    ridge = Ridge(alpha=1.0).fit(train_expression, train_truth)
    test_units = {
        "GSE253719_heldout_donors": {
            "expression": dense(bridge["expression"])[split["test"]],
            "latent": dense(bridge["latent"])[split["test"]],
            "truth": dense(bridge["truth"])[split["test"]],
        },
        "GSE128639_all": {
            "expression": dense(gse128639["expression"]),
            "latent": query_latent["gse128639"],
            "truth": dense(gse128639["protein"])[:, list(gse128639["targets"]).index("HLA_DR_COMPLEX")],
        },
    }
    for condition in np.unique(gse254985["conditions"].astype(str)):
        mask = gse254985["conditions"].astype(str) == condition
        test_units[f"GSE254985_{condition}"] = {
            "expression": dense(gse254985["expression"])[mask],
            "latent": query_latent["gse254985"][mask],
            "truth": dense(gse254985["protein"])[mask, list(gse254985["targets"]).index("HLA_DR_COMPLEX")],
        }

    metric_rows = []
    history_rows = []
    for seed in args.seeds:
        model, history = fine_tune_hla(
            base_model,
            train_expression,
            train_truth,
            train_latent,
            validation_expression,
            validation_truth,
            validation_latent,
            gene_embedding,
            protein_embedding,
            target_index,
            seed=seed,
            epochs=args.epochs,
            patience=args.patience,
            learning_rate=args.learning_rate,
            batch_size=args.batch_size,
        )
        history_rows.extend({"seed": seed, **row} for row in history)
        for unit, block in test_units.items():
            scprotrans_prediction = module.predict_protrans(
                model,
                block["expression"],
                block["latent"],
                gene_embedding,
                protein_embedding,
                batch_size=128,
                device="cpu",
            )[:, target_index]
            predictions = {
                "scprotrans_bridge": scprotrans_prediction,
                "ridge_bridge": ridge.predict(block["expression"]),
                "cognate_rna": module.cognate_rna_score(
                    block["expression"], gene_order=gene_order, target="HLA_DR_COMPLEX"
                ),
            }
            metric_rows.extend(unit_metrics(module, block["truth"], predictions, unit, seed))

    metrics = pd.DataFrame(metric_rows)
    history = pd.DataFrame(history_rows)
    metrics.to_csv(output_dir / "bridge_finetuning_metrics.tsv", sep="\t", index=False)
    history.to_csv(output_dir / "bridge_finetuning_history.tsv", sep="\t", index=False)
    summary = (
        metrics.groupby(["unit", "method"], as_index=False)
        .agg(
            n_cells=("n_cells", "first"),
            mean_spearman=("spearman", "mean"),
            min_spearman=("spearman", "min"),
            max_spearman=("spearman", "max"),
        )
    )
    summary.to_csv(output_dir / "bridge_finetuning_summary.tsv", sep="\t", index=False)

    expanded_reference_latent = np.vstack(
        [reference_latent, dense(bridge["latent"])[split["train"]]]
    )
    expanded_groups = np.concatenate(
        [reference["donors"].astype(str), bridge["donors"].astype(str)[split["train"]]]
    )
    ood_model = module.fit_knn_ood(
        expanded_reference_latent,
        reference_groups=expanded_groups,
        k=10,
        quantile=0.95,
    )
    _, hgsoc_in_domain = module.score_knn_ood(ood_model, query_latent["hgsoc"])
    strict = hgsoc["high_confidence"].astype(bool)
    hgsoc_ood = module.patient_equal_in_domain_fraction(
        hgsoc_in_domain[strict], hgsoc["patients"].astype(str)[strict]
    )

    transfer_decisions = []
    for unit in sorted(test_units):
        unit_summary = summary[summary["unit"] == unit].set_index("method")
        model_metric = float(unit_summary.loc["scprotrans_bridge", "mean_spearman"])
        best_baseline = max(
            float(unit_summary.loc["ridge_bridge", "mean_spearman"]),
            float(unit_summary.loc["cognate_rna", "mean_spearman"]),
        )
        transfer_decisions.append(
            {
                "unit": unit,
                "scprotrans_mean_spearman": model_metric,
                "best_baseline_spearman": best_baseline,
                "delta_vs_best_baseline": model_metric - best_baseline,
                "positive": model_metric > 0,
                "noninferior_margin_0.05": model_metric >= best_baseline - 0.05,
            }
        )
    transfer_passed = all(
        item["positive"] and item["noninferior_margin_0.05"]
        for item in transfer_decisions
    )
    ood_passed = hgsoc_ood["patient_equal_fraction"] >= 0.70
    decision = {
        "route": "exploratory epithelial bridge fine-tuning",
        "split": {
            "train_donors": split["train_donors"].tolist(),
            "validation_donors": split["validation_donors"].tolist(),
            "test_donors": split["test_donors"].tolist(),
        },
        "seeds": args.seeds,
        "transfer_units": transfer_decisions,
        "independent_transfer_passed": transfer_passed,
        "hgsoc_ood_summary": hgsoc_ood,
        "hgsoc_ood_passed_70pct": ood_passed,
        "exploratory_gate_passed": bool(transfer_passed and ood_passed),
        "replaces_original_gate3": False,
        "hgsoc_protein_inference_authorized": False,
        "boundary": (
            "A bridge-trained result is exploratory and cannot replace the frozen preregistered "
            "Gate 3 without a new independent ovarian RNA+ADT validation cohort."
        ),
    }
    (output_dir / "bridge_finetuning_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
