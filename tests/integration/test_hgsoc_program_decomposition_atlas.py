from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "outputs" / "scprotrans_hgsoc_v4" / "tables"


def test_program_atlas_contract_and_hlaii_position():
    atlas = pd.read_csv(TABLE_DIR / "hgsoc_treatment_program_decomposition_atlas.tsv", sep="\t")
    genes = pd.read_csv(TABLE_DIR / "hgsoc_treatment_program_gene_sets.tsv", sep="\t")
    assert atlas.shape[0] == 14
    assert atlas["program"].is_unique
    assert set(atlas["induction_dominance"]) <= {
        "within-dominant",
        "composition-dominant",
        "mixed",
    }
    assert genes.groupby("program")["present"].sum().min() >= 4

    hla = atlas.loc[atlas["program"].eq("HLAII_CD74_CORE")].iloc[0]
    assert hla["estimate_total_change"] > 0
    assert hla["estimate_within_state_component"] > 0
    assert hla["estimate_composition_component"] > 0
    assert hla["bootstrap_positive_fraction_total_change"] == 1
    assert hla["bootstrap_positive_fraction_within_state_component"] == 1
    assert int(hla["absolute_total_rank"]) == 1
    assert (
        atlas.sort_values("estimate_within_state_component", ascending=False)
        .iloc[0]["program"]
        == "HLAII_CD74_CORE"
    )

    ifn = atlas.loc[atlas["program"].eq("IFN_RESPONSE")].iloc[0]
    assert ifn["estimate_total_change"] < 0
    assert ifn["estimate_within_state_component"] < 0
