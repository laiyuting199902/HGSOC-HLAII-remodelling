import csv
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/scprotrans_hgsoc_v4/cycif_hlaii_rescue"


def read_tsv(name: str) -> list[dict[str, str]]:
    with (OUT / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_patient_paired_mhcii_statistics_are_frozen():
    rows = read_tsv("cycif_paired_statistics.tsv")
    lookup = {(row["analysis"], row["metric"]): row for row in rows}

    broad = lookup[("broad_tumor", "mean_mhcii")]
    strict = lookup[("strict_epithelial", "mean_mhcii")]

    assert int(broad["n_patients"]) == 6
    assert int(strict["n_patients"]) == 6
    assert int(broad["positive_n"]) == 5
    assert int(strict["positive_n"]) == 5
    assert float(broad["mean_delta"]) == pytest.approx(0.0838586650, abs=1e-8)
    assert float(strict["mean_delta"]) == pytest.approx(0.0973764863, abs=1e-8)
    assert float(broad["ci_low"]) < 0 < float(broad["ci_high"])
    assert float(strict["ci_low"]) > 0


def test_cross_technology_and_spatial_boundaries_are_frozen():
    cross = {
        row["identity_definition"]: row
        for row in read_tsv("cycif_rna_protein_correlations.tsv")
    }
    assert int(cross["broad_tumor"]["direction_concordant_n"]) == 5
    assert float(cross["broad_tumor"]["spearman_rho"]) == pytest.approx(0.6571428571)

    nearest = read_tsv("cycif_spatial_nearest_neighbor_statistics.tsv")
    assert {row["comparator"] for row in nearest} == {"CD8_T", "myeloid"}
    assert all(int(row["n_patients"]) == 6 for row in nearest)
    assert all(float(row["mean_delta"]) < 0 for row in nearest)
    assert all(float(row["fdr"]) == pytest.approx(0.03125) for row in nearest)

    lineage = {
        row["comparator"]: row
        for row in read_tsv("cycif_lineage_contrast_statistics.tsv")
    }
    for comparator in ("T_cell", "stroma"):
        assert int(lineage[comparator]["positive_n"]) == 6
        assert float(lineage[comparator]["fdr"]) == pytest.approx(0.046875)
