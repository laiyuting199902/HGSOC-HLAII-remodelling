#!/usr/bin/env python3

"""Independent same-site longitudinal validation in GSE201047.

The script deliberately excludes the CD74/HLA-II target genes from epithelial
identity rules. It uses patient-equal, same-site pseudobulk summaries and
reports strict and broad identity definitions as sensitivity analyses.
"""

from __future__ import annotations

import argparse
import json
import re
import tarfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats


CORE_GENES = ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1"]
EPITHELIAL_GENES = [
    "EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "PAX8", "WFDC2",
    "MUC16", "MSLN", "CLDN3", "CLDN4",
]
OVARIAN_GENES = ["PAX8", "WFDC2", "MUC16", "MSLN", "KRT7", "CLDN3", "CLDN4"]
IMMUNE_GENES = ["PTPRC", "CD3D", "CD3E", "MS4A1", "CD79A", "LST1", "TYROBP", "FCER1G", "NKG7"]
STROMAL_GENES = ["COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1", "COL6A2", "C7"]
ENDOTHELIAL_GENES = ["VWF", "KDR", "ENG", "EMCN"]
SELECTED_GENES = sorted(set(CORE_GENES + EPITHELIAL_GENES + IMMUNE_GENES + STROMAL_GENES + ENDOTHELIAL_GENES))

SAMPLE_RE = re.compile(
    r"^(GSM\d+)_(\d+)_([NT])_([A-Z]+)_PT([123])_filtered_gene_bc_matrices_h5\.h5$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        default=str(Path(__file__).resolve().parents[1] / "data" / "raw" / "gse201047"),
    )
    parser.add_argument(
        "--output-dir",
        default=str(Path(__file__).resolve().parents[1] / "outputs" / "scprotrans_hgsoc_v4" / "tables"),
    )
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument(
        "--sensitivity-min-cells",
        type=int,
        nargs="+",
        default=[20, 50, 100],
    )
    parser.add_argument("--bootstrap-iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=260718)
    return parser.parse_args()


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values])


def extract_archive(data_dir: Path) -> list[Path]:
    files = sorted(data_dir.glob("GSM*.h5"))
    if files:
        return files
    archive = data_dir / "GSE201047_RAW.tar"
    if not archive.exists():
        raise FileNotFoundError(f"Missing archive: {archive}")
    with tarfile.open(archive, "r") as handle:
        safe_members = []
        root = data_dir.resolve()
        for member in handle.getmembers():
            target = (data_dir / member.name).resolve()
            if root not in target.parents and target != root:
                raise RuntimeError(f"Unsafe archive member: {member.name}")
            safe_members.append(member)
        handle.extractall(data_dir, members=safe_members, filter="data")
    files = sorted(data_dir.glob("GSM*.h5"))
    if len(files) != 22:
        raise RuntimeError(f"Expected 22 H5 files, found {len(files)}")
    return files


def sample_manifest(files: list[Path]) -> pd.DataFrame:
    rows = []
    for path in files:
        match = SAMPLE_RE.match(path.name)
        if match is None:
            raise ValueError(f"Unexpected sample filename: {path.name}")
        gsm, order, stage, site, patient = match.groups()
        rows.append(
            {
                "gsm": gsm,
                "sample_order": int(order),
                "patient_id": f"PT{patient}",
                "stage": {"N": "chemo-naive", "T": "treated"}[stage],
                "site": site,
                "is_solid_site": site != "A",
                "path": str(path),
            }
        )
    return pd.DataFrame(rows).sort_values("sample_order").reset_index(drop=True)


def read_10x_h5(path: Path) -> tuple[sparse.csc_matrix, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        if "matrix" in handle:
            group = handle["matrix"]
        elif all(key in handle for key in ["data", "indices", "indptr", "shape", "barcodes"]):
            group = handle
        else:
            genome_groups = [
                value
                for value in handle.values()
                if isinstance(value, h5py.Group)
                and all(key in value for key in ["data", "indices", "indptr", "shape", "barcodes"])
            ]
            if len(genome_groups) != 1:
                raise RuntimeError(f"Could not identify a unique 10x matrix group: {path}")
            group = genome_groups[0]
        feature_group = group["features"] if "features" in group else group
        if "name" in feature_group:
            gene_key = "name"
        elif "gene_names" in feature_group:
            gene_key = "gene_names"
        else:
            raise RuntimeError(f"Could not identify gene symbols in 10x H5: {path}")
        genes = decode(feature_group[gene_key][...])
        barcodes = decode(group["barcodes"][...])
        matrix = sparse.csc_matrix(
            (
                group["data"][...],
                group["indices"][...],
                group["indptr"][...],
            ),
            shape=tuple(group["shape"][...]),
        )
    if matrix.shape != (len(genes), len(barcodes)):
        raise RuntimeError(f"H5 dimensions do not match features/barcodes: {path}")
    return matrix, genes, barcodes


def selected_frame(matrix: sparse.csc_matrix, genes: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    lookup = {gene: index for index, gene in enumerate(genes)}
    missing = [gene for gene in SELECTED_GENES if gene not in lookup]
    if missing:
        raise RuntimeError(f"Required genes missing: {', '.join(missing)}")
    indexes = [lookup[gene] for gene in SELECTED_GENES]
    return matrix[indexes, :].toarray().astype(np.float64), {gene: i for i, gene in enumerate(SELECTED_GENES)}


def module_mean(log_expr: np.ndarray, index: dict[str, int], genes: list[str]) -> np.ndarray:
    return log_expr[[index[gene] for gene in genes], :].mean(axis=0)


def detected_count(counts: np.ndarray, index: dict[str, int], genes: list[str]) -> np.ndarray:
    return (counts[[index[gene] for gene in genes], :] > 0).sum(axis=0)


def classify_cells(matrix: sparse.csc_matrix, genes: np.ndarray) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    selected, index = selected_frame(matrix, genes)
    library_size = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    n_features = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(int)
    log_expr = np.log1p(10000.0 * selected / np.maximum(library_size, 1.0))

    epi_score = module_mean(log_expr, index, EPITHELIAL_GENES)
    immune_score = module_mean(log_expr, index, IMMUNE_GENES)
    stromal_score = module_mean(log_expr, index, STROMAL_GENES)
    endothelial_score = module_mean(log_expr, index, ENDOTHELIAL_GENES)
    competitor_score = np.maximum.reduce([immune_score, stromal_score, endothelial_score])
    epi_detect = detected_count(selected, index, EPITHELIAL_GENES)
    ovarian_detect = detected_count(selected, index, OVARIAN_GENES)
    immune_detect = detected_count(selected, index, IMMUNE_GENES)
    ptprc_count = selected[index["PTPRC"], :]

    base_qc = (library_size >= 500) & (n_features >= 200)
    broad = (
        base_qc
        & (epi_detect >= 2)
        & (ovarian_detect >= 1)
        & (epi_score > competitor_score)
        & (ptprc_count == 0)
    )
    strict = (
        base_qc
        & (epi_detect >= 3)
        & (ovarian_detect >= 1)
        & ((epi_score - competitor_score) >= 0.25)
        & (ptprc_count == 0)
        & (immune_detect <= 1)
    )
    audit = pd.DataFrame(
        {
            "library_size": library_size,
            "n_features": n_features,
            "epithelial_score": epi_score,
            "competitor_score": competitor_score,
            "epithelial_margin": epi_score - competitor_score,
            "epithelial_detected": epi_detect,
            "ovarian_detected": ovarian_detect,
            "immune_detected": immune_detect,
            "ptprc_count": ptprc_count,
            "broad_eoc": broad,
            "strict_eoc": strict,
        }
    )
    target_rows = {gene: selected[index[gene], :] for gene in CORE_GENES}
    return audit, target_rows


def pseudobulk_sample(
    matrix: sparse.csc_matrix,
    audit: pd.DataFrame,
    target_rows: dict[str, np.ndarray],
    identity_column: str,
) -> dict[str, float]:
    keep = audit[identity_column].to_numpy(dtype=bool)
    n_cells = int(keep.sum())
    total_library = float(np.asarray(matrix[:, keep].sum()).item()) if n_cells else 0.0
    row: dict[str, float] = {"n_cells": n_cells, "total_library": total_library}
    core_sum = 0.0
    for gene in CORE_GENES:
        gene_sum = float(target_rows[gene][keep].sum()) if n_cells else 0.0
        row[f"count_{gene}"] = gene_sum
        row[f"log2cpm_{gene}"] = np.log2((gene_sum + 0.5) / (total_library + 1.0) * 1e6)
        core_sum += gene_sum
    row["count_core"] = core_sum
    row["log2cpm_core"] = np.log2((core_sum + 0.5 * len(CORE_GENES)) / (total_library + 1.0) * 1e6)
    return row


def patient_same_site_deltas(sample_table: pd.DataFrame, min_cells: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    eligible = sample_table[(sample_table["is_solid_site"]) & (sample_table["n_cells"] >= min_cells)].copy()
    site_counts = (
        eligible.groupby(["identity_definition", "patient_id", "site"])["stage"]
        .nunique()
        .reset_index(name="n_stages")
    )
    common = site_counts[site_counts["n_stages"] == 2][["identity_definition", "patient_id", "site"]]
    matched = eligible.merge(common, on=["identity_definition", "patient_id", "site"], how="inner")
    value_columns = ["log2cpm_core"] + [f"log2cpm_{gene}" for gene in CORE_GENES]
    stage_mean = (
        matched.groupby(["identity_definition", "patient_id", "stage"], as_index=False)[value_columns]
        .mean()
    )
    stage_mean["n_common_sites"] = stage_mean.apply(
        lambda row: int(
            common[
                (common["identity_definition"] == row["identity_definition"])
                & (common["patient_id"] == row["patient_id"])
            ]["site"].nunique()
        ),
        axis=1,
    )
    wide = stage_mean.pivot(index=["identity_definition", "patient_id"], columns="stage", values=value_columns)
    rows = []
    for identity, patient in wide.index:
        row = {"identity_definition": identity, "patient_id": patient}
        sites = common[(common["identity_definition"] == identity) & (common["patient_id"] == patient)]["site"]
        row["common_sites"] = ";".join(sorted(sites.unique()))
        row["n_common_sites"] = int(sites.nunique())
        for value in value_columns:
            pre = float(wide.loc[(identity, patient), (value, "chemo-naive")])
            post = float(wide.loc[(identity, patient), (value, "treated")])
            label = value.replace("log2cpm_", "")
            row[f"pre_{label}"] = pre
            row[f"post_{label}"] = post
            row[f"delta_{label}"] = post - pre
        rows.append(row)
    return matched, pd.DataFrame(rows)


def bootstrap_mean(values: np.ndarray, iterations: int, seed: int) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    draws = rng.choice(values, size=(iterations, len(values)), replace=True).mean(axis=1)
    return float(values.mean()), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def summarize_deltas(deltas: pd.DataFrame, iterations: int, seed: int) -> pd.DataFrame:
    rows = []
    targets = ["core"] + CORE_GENES
    for identity in sorted(deltas["identity_definition"].unique()):
        frame = deltas[deltas["identity_definition"] == identity]
        for offset, target in enumerate(targets):
            values = frame[f"delta_{target}"].to_numpy(dtype=float)
            estimate, low, high = bootstrap_mean(values, iterations, seed + offset)
            positive = int((values > 0).sum())
            nonzero = int((values != 0).sum())
            sign_p = float(stats.binomtest(positive, nonzero, p=0.5, alternative="two-sided").pvalue) if nonzero else np.nan
            rows.append(
                {
                    "identity_definition": identity,
                    "target": target,
                    "n_patients": len(values),
                    "mean_delta_log2cpm": estimate,
                    "bootstrap_ci_low": low,
                    "bootstrap_ci_high": high,
                    "median_delta_log2cpm": float(np.median(values)),
                    "positive_patients": positive,
                    "sign_test_p": sign_p,
                }
            )
    return pd.DataFrame(rows)


def threshold_sensitivity(
    sample_table: pd.DataFrame,
    thresholds: list[int],
    iterations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    delta_frames = []
    summary_frames = []
    for threshold in sorted(set(thresholds)):
        _, deltas = patient_same_site_deltas(sample_table, threshold)
        if deltas.empty:
            continue
        deltas.insert(0, "min_cells", threshold)
        summary = summarize_deltas(deltas, iterations, seed + threshold)
        summary.insert(0, "min_cells", threshold)
        delta_frames.append(deltas)
        summary_frames.append(summary)
    if not delta_frames:
        raise RuntimeError("No threshold sensitivity specification retained a matched patient")
    return pd.concat(delta_frames, ignore_index=True), pd.concat(summary_frames, ignore_index=True)


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    files = extract_archive(data_dir)
    manifest = sample_manifest(files)

    sample_rows = []
    for sample in manifest.itertuples(index=False):
        matrix, genes, barcodes = read_10x_h5(Path(sample.path))
        audit, target_rows = classify_cells(matrix, genes)
        for identity_column, identity_label in [("strict_eoc", "strict"), ("broad_eoc", "broad")]:
            row = sample._asdict()
            row.pop("path")
            row["identity_definition"] = identity_label
            row["n_total_cells"] = matrix.shape[1]
            row.update(pseudobulk_sample(matrix, audit, target_rows, identity_column))
            sample_rows.append(row)

    sample_table = pd.DataFrame(sample_rows)
    matched_sites, deltas = patient_same_site_deltas(sample_table, args.min_cells)
    if deltas.empty:
        raise RuntimeError("No patient had eligible same-site samples at both stages")
    summary = summarize_deltas(deltas, args.bootstrap_iterations, args.seed)
    sensitivity_deltas, sensitivity_summary = threshold_sensitivity(
        sample_table,
        args.sensitivity_min_cells,
        args.bootstrap_iterations,
        args.seed,
    )

    paths = {
        "manifest": output_dir / "gse201047_sample_manifest.tsv",
        "sample_pseudobulk": output_dir / "gse201047_eoc_sample_pseudobulk.tsv",
        "matched_sites": output_dir / "gse201047_same_site_eligible_samples.tsv",
        "patient_deltas": output_dir / "gse201047_same_site_patient_deltas.tsv",
        "summary": output_dir / "gse201047_same_site_summary.tsv",
        "sensitivity_deltas": output_dir / "gse201047_threshold_sensitivity_patient_deltas.tsv",
        "sensitivity_summary": output_dir / "gse201047_threshold_sensitivity_summary.tsv",
        "decision": output_dir / "gse201047_external_validation_decision.json",
    }
    manifest.drop(columns="path").to_csv(paths["manifest"], sep="\t", index=False)
    sample_table.to_csv(paths["sample_pseudobulk"], sep="\t", index=False)
    matched_sites.to_csv(paths["matched_sites"], sep="\t", index=False)
    deltas.to_csv(paths["patient_deltas"], sep="\t", index=False)
    summary.to_csv(paths["summary"], sep="\t", index=False)
    sensitivity_deltas.to_csv(paths["sensitivity_deltas"], sep="\t", index=False)
    sensitivity_summary.to_csv(paths["sensitivity_summary"], sep="\t", index=False)

    strict_core = summary[(summary["identity_definition"] == "strict") & (summary["target"] == "core")].iloc[0]
    broad_core = summary[(summary["identity_definition"] == "broad") & (summary["target"] == "core")].iloc[0]
    decision = {
        "dataset": "GSE201047",
        "analysis": "patient-equal same-solid-site directional replication",
        "min_cells_per_sample": args.min_cells,
        "strict_core_positive_patients": int(strict_core["positive_patients"]),
        "strict_core_n_patients": int(strict_core["n_patients"]),
        "strict_core_mean_delta_log2cpm": float(strict_core["mean_delta_log2cpm"]),
        "broad_core_positive_patients": int(broad_core["positive_patients"]),
        "broad_core_n_patients": int(broad_core["n_patients"]),
        "broad_core_mean_delta_log2cpm": float(broad_core["mean_delta_log2cpm"]),
        "threshold_sensitivity": sensitivity_summary[
            (sensitivity_summary["identity_definition"] == "strict")
            & (sensitivity_summary["target"] == "core")
        ][
            ["min_cells", "n_patients", "mean_delta_log2cpm", "positive_patients"]
        ].to_dict(orient="records"),
        "replication_class": (
            "directionally_concordant"
            if strict_core["positive_patients"] == strict_core["n_patients"]
            and broad_core["positive_patients"] == broad_core["n_patients"]
            else "mixed_or_discordant"
        ),
        "boundary": "Independent but small-n directional replication using target-independent rule-based EOC identity; not a confirmatory validation cohort.",
    }
    paths["decision"].write_text(json.dumps(decision, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(decision, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
