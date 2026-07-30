#!/usr/bin/env python3
"""GSE184880 validation for independent HGSOC protein-layer candidates."""

from __future__ import annotations

import gzip
import importlib.util
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import io, stats


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "geo" / "GSE184880"
DISCOVERY = ROOT / "results" / "scprotrans_hgsoc_independent"
OUT = ROOT / "results" / "scprotrans_hgsoc_independent" / "gse184880_validation"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"


def load_discovery_module():
    script = ROOT / "scripts" / "19_scprotrans_hgsoc_independent_analysis.py"
    spec = importlib.util.spec_from_file_location("scpt_independent", script)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def read_symbols(path: Path) -> list[str]:
    symbols = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            symbols.append(parts[1] if len(parts) > 1 else parts[0])
    return symbols


def sample_prefixes() -> list[str]:
    prefixes = []
    for matrix_file in sorted(RAW.glob("*.matrix.mtx.gz")):
        prefixes.append(matrix_file.name.replace(".matrix.mtx.gz", ""))
    return prefixes


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    mod = load_discovery_module()
    panel = pd.DataFrame(mod.panel_rows())
    genes = sorted(panel["gene"].unique())
    gene_to_module = dict(zip(panel["gene"], panel["module"]))

    sample_rows = []
    gene_rows = []
    pct_rows = []
    for prefix in sample_prefixes():
        status = "Normal" if "_Norm" in prefix else "Cancer"
        genes_file = RAW / f"{prefix}.genes.tsv.gz"
        barcodes_file = RAW / f"{prefix}.barcodes.tsv.gz"
        matrix_file = RAW / f"{prefix}.matrix.mtx.gz"
        symbols = np.array(read_symbols(genes_file), dtype=object)
        with gzip.open(barcodes_file, "rt") as handle:
            n_cells = sum(1 for _ in handle)
        with gzip.open(matrix_file, "rb") as handle:
            mat = io.mmread(handle).tocsc().astype(np.float32)
        total = np.asarray(mat.sum(axis=0)).ravel().astype(np.float32)
        total[total <= 0] = np.nan
        sample_rows.append(
            {
                "sample": prefix,
                "status": status,
                "n_cells": int(n_cells),
                "n_genes": int(mat.shape[0]),
                "matrix_cells": int(mat.shape[1]),
            }
        )
        for gene in genes:
            idx = np.where(symbols == gene)[0]
            if len(idx) == 0:
                gene_rows.append(
                    {
                        "sample": prefix,
                        "status": status,
                        "gene": gene,
                        "module": gene_to_module.get(gene, ""),
                        "present": False,
                        "mean_logcp10k": np.nan,
                    }
                )
                pct_rows.append(
                    {
                        "sample": prefix,
                        "status": status,
                        "gene": gene,
                        "detection_fraction": np.nan,
                    }
                )
                continue
            counts = np.asarray(mat[idx, :].sum(axis=0)).ravel().astype(np.float32)
            norm = np.log1p((counts / total) * 1e4)
            gene_rows.append(
                {
                    "sample": prefix,
                    "status": status,
                    "gene": gene,
                    "module": gene_to_module.get(gene, ""),
                    "present": True,
                    "mean_logcp10k": float(np.nanmean(norm)),
                }
            )
            pct_rows.append(
                {
                    "sample": prefix,
                    "status": status,
                    "gene": gene,
                    "detection_fraction": float(np.nanmean(counts > 0)),
                }
            )

    samples = pd.DataFrame(sample_rows)
    gene_sample = pd.DataFrame(gene_rows)
    pct_sample = pd.DataFrame(pct_rows)
    samples.to_csv(TABLES / "gse184880_sample_summary.tsv", sep="\t", index=False)
    gene_sample.to_csv(TABLES / "gse184880_sample_gene_logcp10k.tsv", sep="\t", index=False)
    pct_sample.to_csv(TABLES / "gse184880_sample_gene_detection.tsv", sep="\t", index=False)

    tests = []
    for gene, sub in gene_sample[gene_sample["present"]].groupby("gene", sort=False):
        normal = sub.loc[sub["status"] == "Normal", "mean_logcp10k"].dropna()
        cancer = sub.loc[sub["status"] == "Cancer", "mean_logcp10k"].dropna()
        if len(normal) >= 2 and len(cancer) >= 2:
            pval = stats.ttest_ind(cancer, normal, equal_var=False).pvalue
            delta = cancer.mean() - normal.mean()
        else:
            pval = np.nan
            delta = np.nan
        tests.append(
            {
                "gene": gene,
                "module": gene_to_module.get(gene, ""),
                "n_normal": len(normal),
                "n_cancer": len(cancer),
                "mean_normal": normal.mean() if len(normal) else np.nan,
                "mean_cancer": cancer.mean() if len(cancer) else np.nan,
                "delta_cancer_minus_normal": delta,
                "welch_p": pval,
            }
        )
    gene_tests = pd.DataFrame(tests)
    gene_tests["welch_fdr"] = mod.bh_adjust(gene_tests["welch_p"])
    gene_tests = gene_tests.sort_values("delta_cancer_minus_normal", ascending=False)
    gene_tests.to_csv(TABLES / "gse184880_gene_cancer_vs_normal.tsv", sep="\t", index=False)

    module_sample_rows = []
    for (sample, status), sub in gene_sample[gene_sample["present"]].groupby(["sample", "status"]):
        row = {"sample": sample, "status": status}
        for module_name, module_sub in sub.groupby("module"):
            row[module_name] = module_sub["mean_logcp10k"].mean()
        module_sample_rows.append(row)
    module_sample = pd.DataFrame(module_sample_rows)
    module_sample.to_csv(TABLES / "gse184880_sample_module_scores.tsv", sep="\t", index=False)

    module_tests = []
    module_cols = [c for c in module_sample.columns if c not in {"sample", "status"}]
    for col in module_cols:
        normal = module_sample.loc[module_sample["status"] == "Normal", col].dropna()
        cancer = module_sample.loc[module_sample["status"] == "Cancer", col].dropna()
        pval = stats.ttest_ind(cancer, normal, equal_var=False).pvalue
        module_tests.append(
            {
                "module": col,
                "n_normal": len(normal),
                "n_cancer": len(cancer),
                "mean_normal": normal.mean(),
                "mean_cancer": cancer.mean(),
                "delta_cancer_minus_normal": cancer.mean() - normal.mean(),
                "welch_p": pval,
            }
        )
    module_tests = pd.DataFrame(module_tests)
    module_tests["welch_fdr"] = mod.bh_adjust(module_tests["welch_p"])
    module_tests = module_tests.sort_values("delta_cancer_minus_normal", ascending=False)
    module_tests.to_csv(TABLES / "gse184880_module_cancer_vs_normal.tsv", sep="\t", index=False)

    # Overlay discovery candidates with validation direction.
    discovery = pd.read_csv(DISCOVERY / "tables" / "ranked_protein_layer_candidates.tsv", sep="\t")
    top_discovery = discovery.head(60)[["cell_type", "feature", "module", "paired_delta_post_minus_naive", "candidate_score"]]
    validation = top_discovery.merge(gene_tests, left_on="feature", right_on="gene", how="left", suffixes=("_discovery", "_validation"))
    validation.to_csv(TABLES / "gse184880_top_discovery_candidate_validation.tsv", sep="\t", index=False)

    # Figures.
    plt.style.use("ggplot")
    fig, ax = plt.subplots(figsize=(8, 5))
    mod_plot = module_tests.sort_values("delta_cancer_minus_normal")
    ax.barh(mod_plot["module"], mod_plot["delta_cancer_minus_normal"], color="#4C78A8")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Cancer minus normal mean logCP10K")
    ax.set_title("GSE184880 module validation")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig6_gse184880_module_validation.png", dpi=300)
    plt.close(fig)

    top_val = validation.dropna(subset=["delta_cancer_minus_normal"]).head(35).copy()
    top_val["label"] = top_val["cell_type"] + ":" + top_val["feature"]
    fig, ax = plt.subplots(figsize=(8.5, max(5, 0.24 * len(top_val))))
    y = np.arange(len(top_val))
    ax.scatter(top_val["delta_cancer_minus_normal"], y, s=70, color="#F58518", label="GSE184880")
    ax.scatter(top_val["paired_delta_post_minus_naive"], y, s=45, color="#54A24B", label="GSE165897 post-NACT")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(top_val["label"])
    ax.invert_yaxis()
    ax.set_xlabel("Delta")
    ax.set_title("Top discovery candidates: chemo delta vs cancer-normal delta")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / "fig7_top_candidate_gse184880_validation.png", dpi=300)
    plt.close(fig)

    summary = OUT / "gse184880_validation_summary.md"
    summary.write_text(
        "\n".join(
            [
                "# GSE184880 validation summary",
                "",
                f"- Samples: {samples.shape[0]}",
                f"- Normal samples: {(samples['status'] == 'Normal').sum()}",
                f"- Cancer samples: {(samples['status'] == 'Cancer').sum()}",
                f"- Tested genes: {gene_tests.shape[0]}",
                "",
                "## Top cancer-up modules",
                "```tsv",
                module_tests.head(10).to_csv(sep="\t", index=False),
                "```",
                "",
                "## Top cancer-up genes",
                "```tsv",
                gene_tests.head(20).to_csv(sep="\t", index=False),
                "```",
                "",
                "## Top discovery candidates with validation",
                "```tsv",
                validation.head(20).to_csv(sep="\t", index=False),
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
