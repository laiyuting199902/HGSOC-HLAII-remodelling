#!/usr/bin/env python3

"""Recover and audit NeoPembrOV sample identities from public source data.

GEO exposes anonymized RNA-seq library identifiers (R200xxx/R210xxx), whereas
the paper source data use clinical subject identifiers.  This script links the
two namespaces with eight expression fingerprints published in Figure 1 and
then verifies the deterministic pre/post library numbering scheme.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata, spearmanr


FINGERPRINTS = [
    "EPCAM",
    "PTPRC (CD45)",
    "CXCL13",
    "JCHAIN",
    "IGHA1-2",
    "IGHG1-4",
    "MZB1 + TNFRSF17",
    "CD8B / FOXP3",
]
COUNT_GENES = {
    "EPCAM",
    "PTPRC",
    "CXCL13",
    "JCHAIN",
    "IGHA1",
    "IGHA2",
    "IGHG1",
    "IGHG2",
    "IGHG3",
    "IGHG4",
    "MZB1",
    "TNFRSF17",
    "CD8B",
    "FOXP3",
}


def parse_args() -> argparse.Namespace:
    root = Path("data/raw/gse227666")
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=root)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project / "outputs" / "scprotrans_hgsoc_v4" / "external_cohort_rescue",
    )
    return parser.parse_args()


def worksheet_frame(path: Path, sheet: str) -> pd.DataFrame:
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = list(workbook[sheet].iter_rows(values_only=True))
    return pd.DataFrame(rows[1:], columns=rows[0])


def source_fingerprints(path: Path) -> pd.DataFrame:
    long = worksheet_frame(path, "Figure 1 a,b")
    long["score"] = pd.to_numeric(long["score"])
    long["subject"] = long["sample_name"].str.replace(r"_(Pre|Post)$", "", regex=True)
    wide = (
        long.pivot_table(
            index=["sample_name", "subject", "Treatment_arm", "time"],
            columns="immune_pop",
            values="score",
        )
        .reset_index()
        .rename_axis(columns=None)
    )
    missing = [name for name in FINGERPRINTS if name not in wide]
    if missing:
        raise RuntimeError(f"Published fingerprints are missing: {missing}")
    wide["CD8B / FOXP3"] = np.log2(wide["CD8B / FOXP3"].clip(lower=1e-6))
    return wide


def geo_metadata(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        text = handle.read()
    rows = []
    pattern = re.compile(r"NACT(\+P)?, (Pre|Post) \((R\d+)\)")
    for section in text.split("^SAMPLE =")[1:]:
        title_match = re.search(r"^!Sample_title = (.*)$", section, re.MULTILINE)
        if title_match is None:
            continue
        match = pattern.search(title_match.group(1))
        if match is None:
            continue
        rows.append(
            {
                "rid": match.group(3),
                "Treatment_arm": "NACT+P" if match.group(1) else "NACT",
                "time": match.group(2),
            }
        )
    metadata = pd.DataFrame(rows)
    if len(metadata) != 109 or metadata["rid"].duplicated().any():
        raise RuntimeError(f"Unexpected GEO metadata dimensions: {metadata.shape}")
    return metadata


def count_fingerprints(path: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    selected: dict[str, np.ndarray] = {}
    with gzip.open(path, "rt") as handle:
        columns = handle.readline().rstrip().split("\t")[1:]
        library_size = np.zeros(len(columns), dtype=float)
        for line in handle:
            fields = line.rstrip().split("\t")
            values = np.asarray(fields[1:], dtype=float)
            library_size += values
            if fields[0] in COUNT_GENES:
                selected[fields[0]] = values
    missing = COUNT_GENES.difference(selected)
    if missing:
        raise RuntimeError(f"Genes needed for record linkage are missing: {sorted(missing)}")

    def log2_cpm(values: np.ndarray) -> np.ndarray:
        return np.log2((values + 0.5) / (library_size + 1.0) * 1e6)

    result = pd.DataFrame({"rid": columns})
    result["EPCAM"] = log2_cpm(selected["EPCAM"])
    result["PTPRC (CD45)"] = log2_cpm(selected["PTPRC"])
    result["CXCL13"] = log2_cpm(selected["CXCL13"])
    result["JCHAIN"] = log2_cpm(selected["JCHAIN"])
    result["IGHA1-2"] = log2_cpm(selected["IGHA1"] + selected["IGHA2"])
    result["IGHG1-4"] = log2_cpm(
        selected["IGHG1"] + selected["IGHG2"] + selected["IGHG3"] + selected["IGHG4"]
    )
    result["MZB1 + TNFRSF17"] = log2_cpm(selected["MZB1"] + selected["TNFRSF17"])
    result["CD8B / FOXP3"] = np.log2(
        (selected["CD8B"] + 0.5) / (selected["FOXP3"] + 0.5)
    )
    return result.merge(metadata, on="rid", validate="one_to_one")


def rank_standardize(frame: pd.DataFrame) -> np.ndarray:
    columns = []
    for feature in FINGERPRINTS:
        values = rankdata(frame[feature].astype(float), method="average")
        standard_deviation = values.std()
        columns.append((values - values.mean()) / max(standard_deviation, 1e-8))
    return np.column_stack(columns)


def assignment_cost(source: pd.DataFrame, counts: pd.DataFrame) -> np.ndarray:
    source_values = rank_standardize(source)
    count_values = rank_standardize(counts)
    return ((source_values[:, None, :] - count_values[None, :, :]) ** 2).mean(axis=2)


def expected_post_id(pre_id: str) -> str:
    number = int(pre_id[1:])
    if pre_id.startswith("R200"):
        return f"R{number + 34}"
    if pre_id.startswith("R210"):
        return f"R{number + 1}"
    raise ValueError(f"Unexpected library identifier: {pre_id}")


def recover_mapping(source: pd.DataFrame, counts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    for arm in ["NACT", "NACT+P"]:
        source_pre = source[(source["Treatment_arm"] == arm) & (source["time"] == "Pre")].reset_index(drop=True)
        count_pre = counts[(counts["Treatment_arm"] == arm) & (counts["time"] == "Pre")].reset_index(drop=True)
        costs = assignment_cost(source_pre, count_pre)
        source_index, count_index = linear_sum_assignment(costs)
        for i, j in zip(source_index, count_index, strict=True):
            ordered = np.sort(costs[i])
            records.append(
                {
                    "subject": source_pre.loc[i, "subject"],
                    "Treatment_arm": arm,
                    "time": "Pre",
                    "rid": count_pre.loc[j, "rid"],
                    "cost": costs[i, j],
                    "next_best_margin": ordered[1] - costs[i, j],
                    "assignment_rank": int(np.argsort(costs[i]).tolist().index(j) + 1),
                    "feature_spearman": spearmanr(
                        source_pre.loc[i, FINGERPRINTS].astype(float),
                        count_pre.loc[j, FINGERPRINTS].astype(float),
                    ).statistic,
                }
            )

    pre_mapping = pd.DataFrame(records)
    post_records = []
    for arm in ["NACT", "NACT+P"]:
        source_post = source[(source["Treatment_arm"] == arm) & (source["time"] == "Post")].reset_index(drop=True)
        count_post = counts[(counts["Treatment_arm"] == arm) & (counts["time"] == "Post")].reset_index(drop=True)
        costs = assignment_cost(source_post, count_post)
        independent_source, independent_count = linear_sum_assignment(costs)
        independent = {
            source_post.loc[i, "subject"]: int(j)
            for i, j in zip(independent_source, independent_count, strict=True)
        }
        arm_pre = pre_mapping[pre_mapping["Treatment_arm"] == arm]
        for row in arm_pre.itertuples(index=False):
            expected = expected_post_id(row.rid)
            source_position = source_post.index[source_post["subject"] == row.subject]
            if len(source_position) != 1 or row.subject not in independent:
                raise RuntimeError(f"Post-treatment source record is unavailable for {row.subject}")
            i, j = int(source_position[0]), independent[row.subject]
            assigned = count_post.loc[j, "rid"]
            ordered = np.sort(costs[i])
            post_records.append(
                {
                    "subject": row.subject,
                    "Treatment_arm": arm,
                    "time": "Post",
                    "rid": assigned,
                    "expected_post_rid": expected,
                    "numbering_consistent": assigned == expected,
                    "cost": costs[i, j],
                    "next_best_margin": ordered[1] - costs[i, j],
                    "assignment_rank": int(np.argsort(costs[i]).tolist().index(j) + 1),
                    "feature_spearman": spearmanr(
                        source_post.loc[i, FINGERPRINTS].astype(float),
                        count_post.loc[j, FINGERPRINTS].astype(float),
                    ).statistic,
                }
            )
    post_mapping = pd.DataFrame(post_records)
    pre_mapping["expected_post_rid"] = pd.NA
    pre_mapping["numbering_consistent"] = pd.NA
    audit = pd.concat([pre_mapping, post_mapping], ignore_index=True)

    pairs = audit.pivot(index=["subject", "Treatment_arm"], columns="time", values="rid").reset_index()
    pairs.columns.name = None
    pairs = pairs.rename(columns={"Pre": "pre_rid", "Post": "post_rid"})
    return audit, pairs


def add_clinical_status(pairs: pd.DataFrame, workbook: Path) -> pd.DataFrame:
    clinical = worksheet_frame(workbook, "Table 1")[
        ["sample_name", "Treatment_arm", "Patient_status"]
    ].drop_duplicates()
    clinical = clinical.rename(columns={"sample_name": "subject"})
    return pairs.merge(clinical, on=["subject", "Treatment_arm"], how="left", validate="one_to_one")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    workbook = args.data_dir / "41467_2024_47000_MOESM4_ESM.xlsx"
    source = source_fingerprints(workbook)
    metadata = geo_metadata(args.data_dir / "GSE227666_family.soft.gz")
    counts = count_fingerprints(args.data_dir / "GSE227666_neopembrov_counts.txt.gz", metadata)
    audit, pairs = recover_mapping(source, counts)
    pairs = add_clinical_status(pairs, workbook)

    post_audit = audit[audit["time"] == "Post"]
    if (
        len(pairs) != 52
        or int(post_audit["numbering_consistent"].sum()) < 51
        or int(audit["assignment_rank"].max()) > 2
        or float(audit["feature_spearman"].median()) < 0.9
    ):
        raise RuntimeError("NeoPembrOV record linkage did not reproduce all 52 public pairs")
    used = set(audit["rid"])
    excluded = metadata[~metadata["rid"].isin(used)].copy()
    excluded["reason"] = "not present in the 52-pair RNA-seq source-data set after published QC"

    audit.to_csv(args.output_dir / "gse227666_record_linkage_audit.tsv", sep="\t", index=False)
    pairs.to_csv(args.output_dir / "gse227666_patient_mapping.tsv", sep="\t", index=False)
    excluded.to_csv(args.output_dir / "gse227666_excluded_libraries.tsv", sep="\t", index=False)
    print(
        f"Recovered {len(pairs)} NeoPembrOV pairs; "
        f"median linkage cost={audit['cost'].median():.4f}; excluded libraries={len(excluded)}"
    )


if __name__ == "__main__":
    main()
