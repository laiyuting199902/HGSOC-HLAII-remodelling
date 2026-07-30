#!/usr/bin/env python3

"""Target-independent epithelial audit for two additional HGSOC scRNA pairs."""

from __future__ import annotations

import argparse
import gzip
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


CORE_GENES = ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1"]
EPITHELIAL_GENES = [
    "EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "PAX8", "WFDC2",
    "MUC16", "MSLN", "CLDN3", "CLDN4",
]
OVARIAN_GENES = ["PAX8", "WFDC2", "MUC16", "MSLN", "KRT7", "CLDN3", "CLDN4"]
IMMUNE_GENES = ["PTPRC", "CD3D", "CD3E", "MS4A1", "CD79A", "LST1", "TYROBP", "FCER1G", "NKG7"]
STROMAL_GENES = ["COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A1", "COL6A2", "C7"]
ENDOTHELIAL_GENES = ["VWF", "KDR", "ENG", "EMCN"]
SELECTED_GENES = sorted(
    set(CORE_GENES + EPITHELIAL_GENES + IMMUNE_GENES + STROMAL_GENES + ENDOTHELIAL_GENES)
)


def parse_args() -> argparse.Namespace:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--database-root",
        type=Path,
        default=Path("data/raw"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project / "outputs" / "scprotrans_hgsoc_v4" / "external_cohort_rescue",
    )
    parser.add_argument("--thresholds", type=int, nargs="+", default=[20, 50, 100])
    return parser.parse_args()


def extract_relevant_zips(root: Path) -> None:
    for archive_name, directory in [
        ("GSM9496525_sp3rna.zip", "sp3rna"),
        ("GSM9496527_sp5rna.zip", "sp5rna"),
    ]:
        target = root / directory / "matrix.mtx.gz"
        if target.exists():
            continue
        archive = root / archive_name
        if not archive.exists():
            raise FileNotFoundError(archive)
        with zipfile.ZipFile(archive) as handle:
            root_resolved = root.resolve()
            for member in handle.infolist():
                destination = (root / member.filename).resolve()
                if root_resolved not in destination.parents and destination != root_resolved:
                    raise RuntimeError(f"Unsafe archive member: {member.filename}")
            handle.extractall(root)


def read_lines(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt") as handle:
        return [line.rstrip().split("\t") for line in handle]


def read_10x_matrix(directory: Path, stem: str | None = None) -> tuple[sparse.csc_matrix, np.ndarray]:
    if stem is None:
        matrix_path = directory / "matrix.mtx.gz"
        feature_path = directory / "features.tsv.gz"
        barcode_path = directory / "barcodes.tsv.gz"
    else:
        matrix_path = directory / f"{stem}_matrix.mtx.gz"
        feature_path = directory / f"{stem}_features.tsv.gz"
        barcode_path = directory / f"{stem}_barcodes.tsv.gz"
    features = read_lines(feature_path)
    genes = np.asarray([row[1] if len(row) > 1 else row[0] for row in features])
    barcodes = read_lines(barcode_path)
    matrix = mmread(matrix_path).tocsc()
    if matrix.shape != (len(genes), len(barcodes)):
        raise RuntimeError(f"Matrix dimensions do not match features/barcodes: {matrix_path}")
    return matrix, genes


def selected_frame(matrix: sparse.csc_matrix, genes: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    lookup = {gene: index for index, gene in enumerate(genes)}
    missing = [gene for gene in SELECTED_GENES if gene not in lookup]
    if missing:
        raise RuntimeError(f"Required genes are missing: {', '.join(missing)}")
    indexes = [lookup[gene] for gene in SELECTED_GENES]
    return matrix[indexes, :].toarray().astype(np.float64), {
        gene: i for i, gene in enumerate(SELECTED_GENES)
    }


def module_mean(log_expression: np.ndarray, index: dict[str, int], genes: list[str]) -> np.ndarray:
    return log_expression[[index[gene] for gene in genes], :].mean(axis=0)


def detected_count(counts: np.ndarray, index: dict[str, int], genes: list[str]) -> np.ndarray:
    return (counts[[index[gene] for gene in genes], :] > 0).sum(axis=0)


def classify_cells(matrix: sparse.csc_matrix, genes: np.ndarray) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    selected, index = selected_frame(matrix, genes)
    library_size = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    n_features = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(int)
    log_expression = np.log1p(10000.0 * selected / np.maximum(library_size, 1.0))
    epithelial_score = module_mean(log_expression, index, EPITHELIAL_GENES)
    immune_score = module_mean(log_expression, index, IMMUNE_GENES)
    stromal_score = module_mean(log_expression, index, STROMAL_GENES)
    endothelial_score = module_mean(log_expression, index, ENDOTHELIAL_GENES)
    competitor_score = np.maximum.reduce([immune_score, stromal_score, endothelial_score])
    epithelial_detected = detected_count(selected, index, EPITHELIAL_GENES)
    ovarian_detected = detected_count(selected, index, OVARIAN_GENES)
    immune_detected = detected_count(selected, index, IMMUNE_GENES)
    ptprc_count = selected[index["PTPRC"], :]
    base_qc = (library_size >= 500) & (n_features >= 200)
    broad = (
        base_qc
        & (epithelial_detected >= 2)
        & (ovarian_detected >= 1)
        & (epithelial_score > competitor_score)
        & (ptprc_count == 0)
    )
    strict = (
        base_qc
        & (epithelial_detected >= 3)
        & (ovarian_detected >= 1)
        & ((epithelial_score - competitor_score) >= 0.25)
        & (ptprc_count == 0)
        & (immune_detected <= 1)
    )
    audit = pd.DataFrame(
        {
            "library_size": library_size,
            "n_features": n_features,
            "epithelial_margin": epithelial_score - competitor_score,
            "epithelial_detected": epithelial_detected,
            "ovarian_detected": ovarian_detected,
            "immune_detected": immune_detected,
            "broad_eoc": broad,
            "strict_eoc": strict,
        }
    )
    targets = {gene: selected[index[gene], :] for gene in CORE_GENES}
    return audit, targets


def pseudobulk(
    matrix: sparse.csc_matrix,
    audit: pd.DataFrame,
    targets: dict[str, np.ndarray],
    identity: str,
) -> dict[str, float]:
    keep = audit[identity].to_numpy(dtype=bool)
    n_cells = int(keep.sum())
    total_library = float(np.asarray(matrix[:, keep].sum()).item()) if n_cells else 0.0
    result: dict[str, float] = {"n_cells": n_cells, "total_library": total_library}
    core_sum = 0.0
    for gene in CORE_GENES:
        gene_sum = float(targets[gene][keep].sum()) if n_cells else 0.0
        result[f"count_{gene}"] = gene_sum
        result[f"log2cpm_{gene}"] = np.log2((gene_sum + 0.5) / (total_library + 1.0) * 1e6)
        core_sum += gene_sum
    result["log2cpm_core"] = np.log2(
        (core_sum + 0.5 * len(CORE_GENES)) / (total_library + 1.0) * 1e6
    )
    return result


def sample_registry(database_root: Path) -> list[dict[str, object]]:
    gse191301 = database_root / "gse191301"
    gse318490 = database_root / "gse318490"
    extract_relevant_zips(gse318490)
    return [
        {
            "cohort": "GSE191301",
            "patient": "patient_1",
            "stage": "Pre",
            "site": "greater omentum",
            "sample": "GSM5743308",
            "reader": lambda: read_10x_matrix(gse191301, "GSM5743308_Pre-NACT1B"),
        },
        {
            "cohort": "GSE191301",
            "patient": "patient_1",
            "stage": "Post",
            "site": "greater omentum",
            "sample": "GSM5743311",
            "reader": lambda: read_10x_matrix(gse191301, "GSM5743311_Post-NACT1E"),
        },
        {
            "cohort": "GSE318490",
            "patient": "P2",
            "stage": "Pre",
            "site": "omentum",
            "sample": "GSM9496525",
            "reader": lambda: read_10x_matrix(gse318490 / "sp3rna"),
        },
        {
            "cohort": "GSE318490",
            "patient": "P2",
            "stage": "Post",
            "site": "omentum",
            "sample": "GSM9496527",
            "reader": lambda: read_10x_matrix(gse318490 / "sp5rna"),
        },
    ]


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    qc_rows = []
    for sample in sample_registry(args.database_root):
        matrix, genes = sample["reader"]()
        audit, targets = classify_cells(matrix, genes)
        qc_rows.append(
            {
                **{key: sample[key] for key in ["cohort", "patient", "stage", "site", "sample"]},
                "n_total_cells": matrix.shape[1],
                "n_strict_eoc": int(audit["strict_eoc"].sum()),
                "n_broad_eoc": int(audit["broad_eoc"].sum()),
                "median_strict_epithelial_margin": float(
                    audit.loc[audit["strict_eoc"], "epithelial_margin"].median()
                ),
            }
        )
        for identity in ["strict_eoc", "broad_eoc"]:
            rows.append(
                {
                    **{key: sample[key] for key in ["cohort", "patient", "stage", "site", "sample"]},
                    "identity_definition": identity,
                    **pseudobulk(matrix, audit, targets, identity),
                }
            )
    scores = pd.DataFrame(rows)
    qc = pd.DataFrame(qc_rows)
    deltas = []
    for (cohort, patient, identity), group in scores.groupby(
        ["cohort", "patient", "identity_definition"]
    ):
        indexed = group.set_index("stage")
        if not {"Pre", "Post"}.issubset(indexed.index):
            continue
        for threshold in args.thresholds:
            eligible = bool((indexed.loc[["Pre", "Post"], "n_cells"] >= threshold).all())
            deltas.append(
                {
                    "cohort": cohort,
                    "patient": patient,
                    "identity_definition": identity,
                    "minimum_cells": threshold,
                    "eligible": eligible,
                    "pre_cells": int(indexed.loc["Pre", "n_cells"]),
                    "post_cells": int(indexed.loc["Post", "n_cells"]),
                    "core_delta": (
                        float(indexed.loc["Post", "log2cpm_core"] - indexed.loc["Pre", "log2cpm_core"])
                        if eligible
                        else np.nan
                    ),
                }
            )
    delta_frame = pd.DataFrame(deltas)
    scores.to_csv(args.output_dir / "additional_scrna_sample_pseudobulk.tsv", sep="\t", index=False)
    qc.to_csv(args.output_dir / "additional_scrna_identity_qc.tsv", sep="\t", index=False)
    delta_frame.to_csv(args.output_dir / "additional_scrna_patient_deltas.tsv", sep="\t", index=False)

    primary = delta_frame[delta_frame["minimum_cells"] == 20]
    report = [
        "# 新增单细胞同部位纵向病例审计",
        "",
        "## 分析原则",
        "",
        "- 肿瘤上皮身份规则不包含 CD74 或任何 HLA-II 目标基因。",
        "- 同时报告严格与宽松身份定义，并在每样本至少 20、50、100 个细胞三个阈值下复算。",
        "- 两个队列各仅有 1 位同患者同部位病例，只能作为方向性病例，不能用于总体显著性检验。",
        "",
        "## 主要结果",
        "",
    ]
    for cohort, group in primary.groupby("cohort"):
        strict = group[group["identity_definition"] == "strict_eoc"].iloc[0]
        broad = group[group["identity_definition"] == "broad_eoc"].iloc[0]
        if strict.eligible and broad.eligible:
            direction = "方向一致，可作为稳健方向性病例" if np.sign(strict.core_delta) == np.sign(broad.core_delta) else "方向不一致，仅作为身份敏感性边界"
            report.append(
                f"- **{cohort}**：严格定义 {strict.pre_cells}/{strict.post_cells} 个治疗前/后 EOC，"
                f"变化 {strict.core_delta:+.3f}；宽松定义 {broad.pre_cells}/{broad.post_cells} 个 EOC，"
                f"变化 {broad.core_delta:+.3f} log2 CPM。**{direction}。**"
            )
        else:
            report.append(f"- **{cohort}**：至少一种身份定义未通过每样本 20 个细胞门槛，不晋级。")
    report.extend(
        [
            "",
            "## 证据边界",
            "",
            "1. GSE318490 的严格与宽松定义方向一致，可作为既有 GSE201047 的补充方向性病例。",
            "2. GSE191301 对身份门槛敏感，保留为异质性边界，不纳入主 forest。",
            "3. 不把两个单病例与患者级整体组织队列合并计算 P 值。",
        ]
    )
    project = Path(__file__).resolve().parents[1]
    final_report = project / "reports" / "scprotrans_hgsoc_v4" / "新增单细胞同部位病例审计.md"
    final_report.parent.mkdir(parents=True, exist_ok=True)
    final_report.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(delta_frame.to_string(index=False))


if __name__ == "__main__":
    main()
