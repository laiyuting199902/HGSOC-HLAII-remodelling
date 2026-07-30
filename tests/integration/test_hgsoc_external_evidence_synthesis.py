import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = ROOT / "outputs/scprotrans_hgsoc_v4/external_cohort_rescue"


def read_tsv(name: str) -> list[dict[str, str]]:
    with (EXTERNAL / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_neopembrov_mapping_and_qc_contract():
    mapping = read_tsv("gse227666_patient_mapping.tsv")
    excluded = read_tsv("gse227666_excluded_libraries.tsv")
    audit = read_tsv("gse227666_record_linkage_audit.tsv")

    assert len(mapping) == 52
    assert len({row["subject"] for row in mapping}) == 52
    assert len({row["pre_rid"] for row in mapping}) == 52
    assert len({row["post_rid"] for row in mapping}) == 52
    assert len(excluded) == 5

    post_audit = [row for row in audit if row["numbering_consistent"]]
    assert len(post_audit) == 52
    assert sum(row["numbering_consistent"].lower() == "true" for row in post_audit) == 51
    assert sum(row["numbering_consistent"].lower() == "false" for row in post_audit) == 1


def test_gse319500_source_split_and_no_double_counting():
    pairs = read_tsv("gse319500_patient_core_deltas.tsv")
    counts: dict[str, int] = {}
    for row in pairs:
        counts[row["source_group"]] = counts.get(row["source_group"], 0) + 1

    assert len(pairs) == 83
    assert counts == {
        "新增 35 对（PMID 32928797）": 35,
        "历史重用：GSE201600（31 对）": 31,
        "历史重用：GSE181597（17 对）": 17,
    }

    evidence = read_tsv("integrated_paired_cohort_evidence.tsv")
    external = [row for row in evidence if row["cohort"] != "GSE266577"]
    assert len(external) == 6
    assert sum(int(row["n_pairs"]) for row in external) == 168
    assert "GSE319500 统一重分析" not in {row["cohort"] for row in evidence}


def test_confirmatory_effects_and_trial_interaction_are_frozen():
    evidence = {row["cohort"]: row for row in read_tsv("integrated_paired_cohort_evidence.tsv")}
    assert math.isclose(float(evidence["GSE201600"]["mean_delta"]), 0.6206802048, rel_tol=1e-8)
    assert float(evidence["GSE201600"]["t_p"]) < 0.01
    assert math.isclose(float(evidence["GSE227666"]["mean_delta"]), 0.6550767884, rel_tol=1e-8)
    assert float(evidence["GSE227666"]["t_p"]) < 0.01

    interaction = read_tsv("gse227666_treatment_interaction.tsv")[0]
    assert abs(float(interaction["estimate"])) < 0.01
    assert float(interaction["p_value"]) > 0.98
