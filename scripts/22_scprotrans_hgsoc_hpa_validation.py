#!/usr/bin/env python3
"""Extract HPA protein-level annotations for HGSOC CD74/HLA-II candidates."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "hpa"
GENE_JSON = RAW / "gene_json"
OUT = ROOT / "results" / "scprotrans_hgsoc_independent" / "hpa_validation"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

GENES = {
    "CD74": "ENSG00000019582",
    "HLA-DRA": "ENSG00000204287",
    "HLA-DRB1": "ENSG00000196126",
    "HLA-DPA1": "ENSG00000231389",
    "HLA-DPB1": "ENSG00000223865",
    "SPP1": "ENSG00000118785",
}


def fetch_missing_json() -> None:
    import subprocess

    GENE_JSON.mkdir(parents=True, exist_ok=True)
    for gene, ens in GENES.items():
        path = GENE_JSON / f"{ens}.json"
        if path.exists() and path.stat().st_size > 0:
            continue
        url = f"https://www.proteinatlas.org/{ens}.json"
        subprocess.run(
            ["curl", "-L", "--fail", "--silent", "--show-error", "-o", str(path), url],
            check=True,
        )


def flatten_top(d: dict | str | None, n: int = 5) -> str:
    if not isinstance(d, dict) or not d:
        return ""
    items = []
    for key, val in d.items():
        try:
            num = float(val)
        except Exception:
            num = float("nan")
        items.append((key, num, val))
    items.sort(key=lambda x: (-1 if pd.isna(x[1]) else -x[1], x[0]))
    return "; ".join(f"{k}:{v}" for k, _, v in items[:n])


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    fetch_missing_json()

    rows = []
    intensity_rows = []
    for gene, ens in GENES.items():
        obj = json.loads((GENE_JSON / f"{ens}.json").read_text())
        ov_tcga = obj.get("Cancer prognostics - Ovary Serous Cystadenocarcinoma (TCGA)", {})
        ov_val = obj.get("Cancer prognostics - Ovary Serous Cystadenocarcinoma (validation)", {})
        rows.append(
            {
                "gene": gene,
                "ensembl": ens,
                "evidence": obj.get("Evidence", ""),
                "hpa_evidence": obj.get("HPA evidence", ""),
                "ih_reliability": obj.get("Reliability (IH)", ""),
                "if_reliability": obj.get("Reliability (IF)", ""),
                "protein_cell_type_specificity": obj.get("Protein cell type specificity", ""),
                "protein_cell_type_distribution": obj.get("Protein cell type distribution", ""),
                "protein_tissue_specificity": obj.get("Protein tissue specificity", ""),
                "protein_tissue_distribution": obj.get("Protein tissue distribution", ""),
                "subcellular_location": ";".join(obj.get("Subcellular location", []) or []),
                "secretome_location": obj.get("Secretome location", ""),
                "secretome_function": obj.get("Secretome function", ""),
                "single_cell_expression_cluster": obj.get("Single cell expression cluster", ""),
                "blood_expression_cluster": obj.get("Blood expression cluster", ""),
                "tissue_expression_cluster": obj.get("Tissue expression cluster", ""),
                "rna_cancer_specificity": obj.get("RNA cancer specificity", ""),
                "rna_cancer_distribution": obj.get("RNA cancer distribution", ""),
                "ov_tcga_prognostic": ov_tcga.get("prognostic", ""),
                "ov_tcga_is_prognostic": ov_tcga.get("is_prognostic", ""),
                "ov_tcga_p": ov_tcga.get("p_val", ""),
                "ov_validation_prognostic": ov_val.get("prognostic", ""),
                "ov_validation_is_prognostic": ov_val.get("is_prognostic", ""),
                "ov_validation_p": ov_val.get("p_val", ""),
                "top_protein_cell_type_intensity": flatten_top(obj.get("Protein cell type specific Intensity"), 5),
                "top_protein_tissue_intensity": flatten_top(obj.get("Protein tissue specific Intensity"), 5),
                "top_rna_single_cell_group": flatten_top(obj.get("RNA single cell type group specific nCPM"), 5),
            }
        )
        intens = obj.get("Protein cell type specific Intensity")
        if isinstance(intens, dict):
            for cell_type, value in intens.items():
                try:
                    numeric = float(value)
                except Exception:
                    numeric = None
                intensity_rows.append(
                    {
                        "gene": gene,
                        "cell_type": cell_type,
                        "protein_cell_type_intensity": numeric,
                    }
                )

    summary = pd.DataFrame(rows)
    summary.to_csv(TABLES / "hpa_candidate_summary.tsv", sep="\t", index=False)
    intensity = pd.DataFrame(intensity_rows)
    intensity.to_csv(TABLES / "hpa_protein_cell_type_intensity.tsv", sep="\t", index=False)

    # Plot top protein cell type intensities for genes with available protein cell-type data.
    if not intensity.empty:
        top = intensity.sort_values("protein_cell_type_intensity", ascending=False).groupby("gene").head(5)
        top["label"] = top["gene"] + ":" + top["cell_type"]
        fig, ax = plt.subplots(figsize=(9, max(4.8, 0.28 * len(top))))
        y = range(len(top))
        ax.barh(list(y), top["protein_cell_type_intensity"], color="#4C78A8")
        ax.set_yticks(list(y))
        ax.set_yticklabels(top["label"])
        ax.invert_yaxis()
        ax.set_xlabel("HPA protein cell-type intensity")
        ax.set_title("HPA protein-level cell-type annotations")
        fig.tight_layout()
        fig.savefig(FIGURES / "hpa_candidate_protein_cell_type_intensity.png", dpi=300)
        plt.close(fig)

    report = OUT / "hpa_validation_summary.md"
    report.write_text(
        "\n".join(
            [
                "# HPA candidate validation summary",
                "",
                "Source: Human Protein Atlas `proteinatlas.tsv.zip` and per-gene JSON endpoints.",
                "",
                "## Candidate summary",
                "```tsv",
                summary.to_csv(sep="\t", index=False),
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
