#!/usr/bin/env python3

"""Frozen scProTrans transfer audit in GSE253719 epithelial RNA+ADT cells."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "39_scprotrans_crossdomain_training.py"


def load_training_module():
    spec = importlib.util.spec_from_file_location("scprotrans_crossdomain", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {MODULE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--h5ad",
        default=(
            "data/raw/gse253719/"
            "GSE253719_Biopsy_RNAADT_Final.h5ad"
        ),
    )
    parser.add_argument(
        "--application-dir",
        default=str(
            ROOT
            / "outputs"
            / "scprotrans_hgsoc_v4"
            / "scprotrans_benchmark"
            / "pilot"
            / "h500_sv10_pt10_seed260716"
            / "application"
        ),
    )
    parser.add_argument(
        "--prepared-dir",
        default=str(
            ROOT
            / "outputs"
            / "scprotrans_hgsoc_v4"
            / "scprotrans_benchmark"
            / "prepared"
            / "g128_r250_h2600_s260716"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=str(
            ROOT
            / "outputs"
            / "scprotrans_hgsoc_v4"
            / "scprotrans_benchmark"
            / "gse253719_bridge"
        ),
    )
    parser.add_argument("--max-cells-per-donor", type=int, default=250)
    parser.add_argument("--seed", type=int, default=260718)
    return parser.parse_args()


def sample_epithelial_cells(obs: pd.DataFrame, maximum: int, seed: int) -> np.ndarray:
    epithelial = obs["220414 JH COARSE"].astype(str).to_numpy() == "Epithelial"
    donors = obs["CoLabs_patient"].astype(str).to_numpy()
    rng = np.random.default_rng(seed)
    retained = []
    for donor in sorted(np.unique(donors[epithelial])):
        candidates = np.flatnonzero(epithelial & (donors == donor))
        if candidates.size > maximum:
            candidates = rng.choice(candidates, size=maximum, replace=False)
        retained.extend(np.asarray(candidates, dtype=int).tolist())
    return np.asarray(sorted(retained), dtype=int)


def donor_metrics(module, truth, predictions, donors) -> pd.DataFrame:
    rows = []
    for method, values in predictions.items():
        for donor in sorted(np.unique(donors)):
            mask = donors == donor
            metrics = module.compute_prediction_metrics(truth[mask], values[mask])
            rows.append({"method": method, "donor": donor, "n_cells": int(mask.sum()), **metrics})
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    module = load_training_module()
    h5ad_path = Path(args.h5ad)
    application_dir = Path(args.application_dir)
    prepared_dir = Path(args.prepared_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model, metadata = module.load_model_bundle(application_dir / "scprotrans_model_bundle.pt")
    gene_order = list(metadata["gene_order"])
    protein_order = list(metadata["protein_order"])
    target_index = protein_order.index("HLA_DR_COMPLEX")

    backed = ad.read_h5ad(h5ad_path, backed="r")
    try:
        missing = [gene for gene in gene_order if gene not in backed.var_names]
        if "HLA-DR" not in backed.var_names or backed.var.loc["HLA-DR", "assay"] != "ADT":
            raise RuntimeError("Measured HLA-DR ADT was not found")
        rows = sample_epithelial_cells(backed.obs, args.max_cells_per_donor, args.seed)
        present_genes = [gene for gene in gene_order if gene in backed.var_names]
        columns = [backed.var_names.get_loc(gene) for gene in present_genes] + [
            backed.var_names.get_loc("HLA-DR")
        ]
        subset = backed[rows, columns].to_memory()
        counts = sparse.csr_matrix(subset.layers["counts"])
        present_counts = counts[:, : len(present_genes)].tocsr()
        present_lookup = {gene: index for index, gene in enumerate(present_genes)}
        zero_column = sparse.csr_matrix((present_counts.shape[0], 1), dtype=present_counts.dtype)
        rna_counts = sparse.hstack(
            [
                present_counts[:, present_lookup[gene]]
                if gene in present_lookup
                else zero_column
                for gene in gene_order
            ],
            format="csr",
        )
        truth = module.transform_target_adt_counts(
            counts[:, len(present_genes)].toarray().reshape(-1, 1)
        ).reshape(-1)
        donors = subset.obs["CoLabs_patient"].astype(str).to_numpy()
        diseases = subset.obs["disease"].astype(str).to_numpy()
        batches = subset.obs["LIBRARY"].astype(str).to_numpy()
        cells = subset.obs_names.astype(str).to_numpy()
    finally:
        backed.file.close()

    expression = module.normalize_total_log1p(rna_counts)
    latent = module.map_scvi_query_latent(
        rna_counts,
        batches=batches,
        gene_order=gene_order,
        reference_model_dir=application_dir / "scvi_reference",
        output_model_dir=output_dir / "scvi_gse253719",
        max_epochs=0,
        seed=args.seed,
        progress=False,
    )
    gene_embedding = np.load(prepared_dir / "gene_embedding.npy")
    protein_embedding = np.load(prepared_dir / "protein_embedding.npy")
    scprotrans = module.predict_protrans(
        model,
        expression,
        latent,
        gene_embedding,
        protein_embedding,
        batch_size=128,
        device="cpu",
    )[:, target_index]

    reference = module.load_prepared_block(prepared_dir, "reference")
    ridge = module.fit_ridge_baseline(reference["expression"], reference["protein"], alpha=1.0)
    ridge_prediction = module.predict_ridge_baseline(ridge, expression)[:, target_index]
    cognate = module.cognate_rna_score(
        expression, gene_order=gene_order, target="HLA_DR_COMPLEX"
    )
    predictions = {
        "scprotrans": scprotrans,
        "ridge": ridge_prediction,
        "cognate_rna": cognate,
    }

    pooled_rows = []
    for method, values in predictions.items():
        pooled_rows.append(
            {
                "method": method,
                "n_cells": len(truth),
                "n_donors": int(np.unique(donors).size),
                **module.compute_prediction_metrics(truth, values),
            }
        )
    pooled = pd.DataFrame(pooled_rows)
    per_donor = donor_metrics(module, truth, predictions, donors)
    donor_summary = (
        per_donor.groupby("method", as_index=False)
        .agg(
            n_donors=("donor", "nunique"),
            finite_donors=("spearman", lambda x: int(np.isfinite(x).sum())),
            donor_equal_spearman=("spearman", "mean"),
        )
    )

    with np.load(application_dir / "knn_ood_model.npz") as archive:
        ood_model = {
            "reference_latent": archive["reference_latent"],
            "k": int(archive["k"][0]),
            "quantile": float(archive["quantile"][0]),
            "threshold": float(archive["threshold"][0]),
        }
    ood_score, in_domain = module.score_knn_ood(ood_model, latent)
    donor_domain = pd.DataFrame({"donor": donors, "in_domain": in_domain}).groupby(
        "donor", as_index=False
    )["in_domain"].mean()

    cell_table = pd.DataFrame(
        {
            "cell": cells,
            "donor": donors,
            "disease": diseases,
            "hla_dr_adt_log1p": truth,
            "scprotrans": scprotrans,
            "ridge": ridge_prediction,
            "cognate_rna": cognate,
            "ood_score": ood_score,
            "in_domain": in_domain,
        }
    )
    pooled.to_csv(output_dir / "gse253719_pooled_metrics.tsv", sep="\t", index=False)
    per_donor.to_csv(output_dir / "gse253719_donor_metrics.tsv", sep="\t", index=False)
    donor_summary.to_csv(output_dir / "gse253719_donor_summary.tsv", sep="\t", index=False)
    donor_domain.to_csv(output_dir / "gse253719_donor_ood.tsv", sep="\t", index=False)
    cell_table.to_csv(
        output_dir / "gse253719_epithelial_predictions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    np.savez_compressed(
        output_dir / "gse253719_bridge_arrays.npz",
        expression=np.asarray(expression.toarray() if sparse.issparse(expression) else expression, dtype=np.float32),
        latent=np.asarray(latent, dtype=np.float32),
        truth=np.asarray(truth, dtype=np.float32),
        donors=np.asarray(donors, dtype=str),
        diseases=np.asarray(diseases, dtype=str),
        cells=np.asarray(cells, dtype=str),
        in_domain=np.asarray(in_domain, dtype=bool),
    )

    pooled_lookup = pooled.set_index("method")["spearman"].to_dict()
    donor_lookup = donor_summary.set_index("method")["donor_equal_spearman"].to_dict()
    pooled_best_baseline = max(pooled_lookup["ridge"], pooled_lookup["cognate_rna"])
    donor_best_baseline = max(donor_lookup["ridge"], donor_lookup["cognate_rna"])
    decision = {
        "dataset": "GSE253719",
        "role": "frozen independent epithelial bridge test",
        "cell_count": int(len(truth)),
        "donor_count": int(np.unique(donors).size),
        "missing_frozen_genes_zero_filled": missing,
        "scprotrans_pooled_spearman": float(pooled_lookup["scprotrans"]),
        "best_baseline_pooled_spearman": float(pooled_best_baseline),
        "pooled_delta_vs_best_baseline": float(
            pooled_lookup["scprotrans"] - pooled_best_baseline
        ),
        "scprotrans_donor_equal_spearman": float(donor_lookup["scprotrans"]),
        "best_baseline_donor_equal_spearman": float(donor_best_baseline),
        "donor_equal_delta_vs_best_baseline": float(
            donor_lookup["scprotrans"] - donor_best_baseline
        ),
        "cell_in_domain_fraction": float(np.mean(in_domain)),
        "donor_equal_in_domain_fraction": float(donor_domain["in_domain"].mean()),
        "noninferior_pooled_margin_0.05": bool(
            pooled_lookup["scprotrans"] >= pooled_best_baseline - 0.05
        ),
        "noninferior_donor_equal_margin_0.05": bool(
            donor_lookup["scprotrans"] >= donor_best_baseline - 0.05
        ),
        "boundary": (
            "This non-ovarian epithelial RNA+ADT dataset tests frozen cross-epithelial "
            "transfer only and cannot validate longitudinal HGSOC protein change."
        ),
    }
    (output_dir / "gse253719_bridge_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
