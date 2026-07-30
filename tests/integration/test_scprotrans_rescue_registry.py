"""Regression tests for the frozen scProTrans rescue dataset registry."""

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REGISTRY = PROJECT_ROOT / "outputs/scprotrans_hgsoc_v4/tables/scprotrans_rescue_dataset_registry.tsv"
REPORT = PROJECT_ROOT / "reports/scprotrans_hgsoc_rescue_data_audit_2026-07-17.md"

REQUIRED_COLUMNS = {
    "branch",
    "resource_name",
    "resource_type",
    "accession_or_id",
    "disease_context",
    "rna_protein_capable",
    "accessible_now",
    "failure_reason",
    "promotion_status",
}
ALLOWED_PROMOTION_STATUSES = {"candidate", "blocked", "not_found", "passed"}
REQUIRED_AUDIT_COLUMNS = {
    "resource_discovery_status",
    "access_status",
    "modality_scope",
    "rescue_eligibility",
    "checked_at",
    "source_url_or_query",
    "evidence_note",
}
ALLOWED_DISCOVERY_STATUSES = {"identified", "not_found"}
ALLOWED_ACCESS_STATUSES = {"accessible", "restricted", "not_applicable"}
ALLOWED_MODALITY_SCOPES = {
    "not_available",
    "measured_spatial_protein",
    "single_cell_rna_only",
    "bulk_matched_rna_protein_context",
}
ALLOWED_RESCUE_ELIGIBILITIES = {"ineligible", "pending_access", "eligible"}
KEY_TEXT_COLUMNS = {
    "branch",
    "resource_name",
    "resource_type",
    "accession_or_id",
    "disease_context",
    "failure_reason",
    "promotion_status",
    "resource_discovery_status",
    "access_status",
    "modality_scope",
    "rescue_eligibility",
    "checked_at",
    "source_url_or_query",
    "evidence_note",
}


def read_registry() -> pd.DataFrame:
    return pd.read_csv(REGISTRY, sep="\t")


def test_registry_has_required_columns_status_domains_and_audit_trail():
    df = read_registry()

    assert REQUIRED_COLUMNS.issubset(df.columns)
    assert REQUIRED_AUDIT_COLUMNS.issubset(df.columns)
    assert len(df) == 4, "The frozen rescue registry must contain exactly four records."
    assert df["accession_or_id"].is_unique
    for column in KEY_TEXT_COLUMNS:
        assert df[column].notna().all()
        assert (df[column].astype(str).str.strip() != "").all()

    assert set(df["promotion_status"]).issubset(ALLOWED_PROMOTION_STATUSES)
    assert set(df["resource_discovery_status"]).issubset(ALLOWED_DISCOVERY_STATUSES)
    assert set(df["access_status"]).issubset(ALLOWED_ACCESS_STATUSES)
    assert set(df["modality_scope"]).issubset(ALLOWED_MODALITY_SCOPES)
    assert set(df["rescue_eligibility"]).issubset(ALLOWED_RESCUE_ELIGIBILITIES)
    assert df["rna_protein_capable"].isin([True, False]).all()
    assert df["accessible_now"].isin([True, False]).all()


def test_registry_contains_current_known_rescue_branches_without_promotion():
    df = read_registry().set_index("accession_or_id", drop=False)

    ovarian_search = df.loc["none_found"]
    assert ovarian_search["branch"] == "A"
    assert ovarian_search["disease_context"] == "human ovarian/HGSOC"
    assert ovarian_search["rna_protein_capable"] == True
    assert ovarian_search["accessible_now"] == False
    assert ovarian_search["promotion_status"] == "blocked"
    assert ovarian_search["resource_discovery_status"] == "not_found"
    assert ovarian_search["access_status"] == "not_applicable"
    assert ovarian_search["modality_scope"] == "not_available"
    assert ovarian_search["rescue_eligibility"] == "ineligible"
    assert ovarian_search["checked_at"] == "2026-07-17"
    assert "ncbi.nlm.nih.gov" in ovarian_search["source_url_or_query"]
    assert "no clean public hit" in ovarian_search["evidence_note"]

    spatial_atlas = df.loc["syn72380119"]
    assert spatial_atlas["branch"] == "B"
    assert spatial_atlas["accessible_now"] == False
    assert "region restriction" in spatial_atlas["failure_reason"]
    assert spatial_atlas["promotion_status"] == "candidate"
    assert spatial_atlas["resource_discovery_status"] == "identified"
    assert spatial_atlas["access_status"] == "restricted"
    assert spatial_atlas["modality_scope"] == "measured_spatial_protein"
    assert spatial_atlas["rescue_eligibility"] == "pending_access"

    gse184880 = df.loc["GSE184880"]
    assert gse184880["accessible_now"] == True
    assert gse184880["rna_protein_capable"] == False
    assert "no protein modality" in gse184880["failure_reason"]
    assert gse184880["promotion_status"] == "not_found"
    assert gse184880["resource_discovery_status"] == "identified"
    assert gse184880["access_status"] == "accessible"
    assert gse184880["modality_scope"] == "single_cell_rna_only"
    assert gse184880["rescue_eligibility"] == "ineligible"

    bulk_context = df.loc["TCGA-OV+CPTAC"]
    assert bulk_context["accessible_now"] == True
    assert bulk_context["rna_protein_capable"] == True
    assert "not tumour-cell longitudinal protein" in bulk_context["failure_reason"]
    assert bulk_context["promotion_status"] == "blocked"
    assert bulk_context["resource_discovery_status"] == "identified"
    assert bulk_context["access_status"] == "accessible"
    assert bulk_context["modality_scope"] == "bulk_matched_rna_protein_context"
    assert bulk_context["rescue_eligibility"] == "ineligible"

    assert not (df["promotion_status"] == "passed").any()


def test_audit_report_locks_date_mainline_boundary_and_candidate_language():
    text = REPORT.read_text(encoding="utf-8")

    assert "2026-07-17" in text
    assert "`scProTrans` 仍不能升级为主结论" in text
    assert "不得将 `candidate` 表述为已通过 rescue" in text
