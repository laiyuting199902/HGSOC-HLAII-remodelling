from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
EXTERNAL = ROOT / "outputs" / "scprotrans_hgsoc_v4" / "external_cohort_rescue"
REPORT = ROOT / "reports" / "scprotrans_hgsoc_v4" / "GSE143897患者配对外部纵向审计.md"


def test_gse143897_recovers_the_publicly_linked_18_patient_pairs():
    manifest = pd.read_csv(EXTERNAL / "gse143897_pair_manifest.tsv", sep="\t")
    assert manifest.shape[0] == 36
    assert manifest["patient_id"].nunique() == 18
    assert manifest.groupby("patient_id")["stage"].agg(set).eq({"Pre", "Post"}).all()
    assert manifest["sample"].is_unique


def test_gse143897_core_effect_is_patient_level_and_bulk_bounded():
    summary = pd.read_csv(EXTERNAL / "gse143897_longitudinal_summary.tsv", sep="\t")
    row = summary.iloc[0]
    assert row["cohort"] == "GSE143897"
    assert row["n_pairs"] == 18
    assert row["positive_pairs"] == 14
    assert row["mean_delta"] > 0
    assert row["t_p"] < 0.05

    report = REPORT.read_text(encoding="utf-8")
    assert "整体组织 RNA" in report
    assert "不能将结果定位为肿瘤细胞内在变化" in report
    assert "不能替代单细胞恶性身份或 DNA 基因型锚定" in report


def test_gse143897_program_matrix_keeps_core_coverage_and_ifn_boundary():
    programs = pd.read_csv(EXTERNAL / "gse143897_program_summary.tsv", sep="\t")
    assert programs.shape[0] == 14
    core = programs.loc[programs["program"].eq("HLAII_CD74_CORE")].iloc[0]
    assert core["n_genes_requested"] == 5
    assert core["n_genes_detected"] == 5
    assert core["mean_delta"] > 0
    ifn = programs.loc[programs["program"].eq("IFN_RESPONSE")].iloc[0]
    assert ifn["mean_delta"] < core["mean_delta"]
