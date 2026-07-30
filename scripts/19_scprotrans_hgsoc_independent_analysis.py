#!/usr/bin/env python3
"""Independent HGSOC scProTrans-ready protein-layer analysis.

This script uses public GSE165897 HGSOC scRNA-seq data as a discovery
dataset and produces patient-paired treatment-naive vs post-NACT summaries
for a biologically curated protein-coding panel. It also audits the local
scProTrans installation so the report can distinguish formal scProTrans
requirements from this RNA-derived protein-layer baseline.
"""

from __future__ import annotations

import gzip
import math
import os
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
CELL_INFO = Path(
    "data/raw/gse165897/"
    "GSE165897_cellInfo_HGSOC.tsv.gz"
)
COUNTS = Path(
    "data/raw/gse165897/"
    "GSE165897_UMIcounts_HGSOC.tsv.gz"
)
PROTRANS_DIR = ROOT / "external" / "ProTrans"
OUT = ROOT / "results" / "scprotrans_hgsoc_independent"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(p)
    out = np.full_like(p, np.nan, dtype=float)
    if ok.sum() == 0:
        return pd.Series(out, index=pvalues.index)
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    restored = np.empty_like(adj)
    restored[order] = adj
    out[ok] = restored
    return pd.Series(out, index=pvalues.index)


def panel_rows() -> list[dict[str, str]]:
    modules = {
        "Antigen_presentation": [
            "B2M",
            "HLA-A",
            "HLA-B",
            "HLA-C",
            "HLA-DRA",
            "HLA-DRB1",
            "HLA-DPA1",
            "HLA-DPB1",
            "TAP1",
            "TAP2",
            "CD74",
        ],
        "Immune_checkpoint": [
            "CD274",
            "PDCD1",
            "CTLA4",
            "LAG3",
            "TIGIT",
            "HAVCR2",
            "CD47",
            "SIRPA",
            "LGALS9",
            "VSIR",
        ],
        "Cytokine_JAK_STAT": [
            "IL6",
            "IL6ST",
            "JAK1",
            "JAK2",
            "STAT1",
            "STAT3",
            "TNF",
            "TNFRSF1A",
            "TNFRSF1B",
            "CXCL8",
            "CXCR2",
            "CCL2",
            "CCR2",
            "CCL5",
            "CCR5",
            "CXCL12",
            "CXCR4",
            "CXCL16",
            "CXCR6",
        ],
        "TGF_beta_VEGF": [
            "TGFB1",
            "TGFB2",
            "TGFBR1",
            "TGFBR2",
            "VEGFA",
            "VEGFB",
            "VEGFC",
            "KDR",
            "FLT1",
            "ENG",
            "ANGPT2",
        ],
        "ECM_adhesion": [
            "COL1A1",
            "COL1A2",
            "COL3A1",
            "COL5A1",
            "COL5A2",
            "COL6A1",
            "COL6A2",
            "FN1",
            "THBS1",
            "SPARC",
            "MMP2",
            "MMP9",
            "MMP14",
            "TIMP1",
            "ITGA5",
            "ITGAV",
            "ITGB1",
            "ITGB2",
            "ICAM1",
            "VCAM1",
            "CD44",
        ],
        "Epithelial_tumor_surface": [
            "EPCAM",
            "MSLN",
            "MUC16",
            "CLDN3",
            "CLDN4",
            "KRT8",
            "KRT18",
            "KRT19",
            "CDH1",
            "CDH2",
            "VIM",
            "PAX8",
        ],
        "CAF_pericyte_endothelial": [
            "DCN",
            "LUM",
            "FAP",
            "ACTA2",
            "TAGLN",
            "PDGFRA",
            "PDGFRB",
            "RGS5",
            "MYH11",
            "PECAM1",
            "VWF",
            "ESAM",
        ],
        "Myeloid_macrophage": [
            "LYZ",
            "LST1",
            "CD68",
            "CD14",
            "FCGR3A",
            "MS4A7",
            "CSF1",
            "CSF1R",
            "SPP1",
            "MIF",
            "NLRP3",
        ],
        "T_NK_cytotoxic": [
            "CD3D",
            "CD3E",
            "CD4",
            "CD8A",
            "CD8B",
            "NKG7",
            "GNLY",
            "PRF1",
            "GZMB",
            "GZMA",
            "IFNG",
        ],
        "Stress_transport_translation": [
            "SLC7A1",
            "SLC7A5",
            "SLC2A1",
            "SLC9A1",
            "SLC16A1",
            "ATF4",
            "DDIT3",
            "HSPA5",
            "HSP90B1",
            "XBP1",
            "MTOR",
            "EIF4E",
            "RPS6",
            "RPLP0",
        ],
        "Mitochondrial_OXPHOS": [
            "MDH2",
            "SDHA",
            "IDH3A",
            "NDUFA9",
            "NDUFS1",
            "UQCRC1",
            "COX4I1",
            "COX5A",
            "ATP5F1A",
            "ATP5MC1",
            "TOMM20",
        ],
        "DNA_damage_apoptosis": [
            "PARP1",
            "PARP4",
            "BCL2",
            "BAX",
            "MCL1",
            "CASP3",
            "CASP8",
            "ERCC1",
            "BRCA1",
            "BRCA2",
            "RAD51",
        ],
    }
    rows = []
    seen = set()
    for module, genes in modules.items():
        for gene in genes:
            if gene not in seen:
                seen.add(gene)
                rows.append({"gene": gene, "module": module})
    return rows


def lr_pairs() -> pd.DataFrame:
    rows = [
        ("CXCL12", "CXCR4", "CAF/tumor chemotaxis"),
        ("CXCL16", "CXCR6", "T-cell recruitment"),
        ("SPP1", "CD44", "myeloid/ECM tumor support"),
        ("MIF", "CD74", "myeloid inflammatory signaling"),
        ("CD47", "SIRPA", "anti-phagocytosis"),
        ("LGALS9", "HAVCR2", "T-cell exhaustion"),
        ("CD274", "PDCD1", "PD-1/PD-L1 checkpoint"),
        ("TGFB1", "TGFBR2", "TGF-beta signaling"),
        ("VEGFA", "KDR", "angiogenesis"),
        ("THBS1", "CD47", "matrix immune checkpoint"),
        ("FN1", "ITGA5", "fibronectin-integrin"),
        ("COL1A1", "ITGB1", "collagen-integrin"),
        ("ICAM1", "ITGB2", "adhesion immune recruitment"),
        ("CCL2", "CCR2", "monocyte recruitment"),
        ("CCL5", "CCR5", "T/myeloid chemotaxis"),
        ("IL6", "IL6ST", "IL6/JAK-STAT"),
        ("TNF", "TNFRSF1B", "TNF receptor signaling"),
        ("CSF1", "CSF1R", "macrophage survival"),
    ]
    return pd.DataFrame(rows, columns=["ligand", "receptor", "axis"])


def audit_scprotrans() -> pd.DataFrame:
    rows = []
    paths = {
        "repo": PROTRANS_DIR,
        "main_script": PROTRANS_DIR / "code" / "ProTrans.py",
        "gene_embedding": PROTRANS_DIR / "dataset" / "gene_embedding" / "dna2vec_1w.npz",
        "protein_embedding": PROTRANS_DIR / "dataset" / "protein_embedding" / "embedding_ProtT5.h5",
        "protein_embedding_alt": PROTRANS_DIR / "dataset" / "protein_embedding" / "per-protein.h5",
    }
    for name, path in paths.items():
        rows.append(
            {
                "item": name,
                "path": str(path),
                "exists": path.exists(),
                "note": "",
            }
        )
    try:
        rev = subprocess.check_output(
            ["git", "-C", str(PROTRANS_DIR), "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception as exc:  # pragma: no cover - audit only
        rev = f"NA: {exc}"
    rows.append({"item": "git_revision", "path": rev, "exists": True, "note": ""})
    rows.append(
        {
            "item": "formal_prediction_status",
            "path": "",
            "exists": False,
            "note": (
                "Official ProTrans.py requires paired rna.csv and protein.csv to train/evaluate; "
                "the cloned repository does not include a pretrained model for direct HGSOC scRNA inference."
            ),
        }
    )
    return pd.DataFrame(rows)


def read_selected_counts(counts_path: Path, selected_genes: list[str]) -> tuple[pd.DataFrame, list[str]]:
    selected = set(selected_genes)
    rows = {}
    cells: list[str] | None = None
    with gzip.open(counts_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        for line in handle:
            if not line:
                continue
            first_tab = line.find("\t")
            gene = line[:first_tab]
            if gene not in selected:
                continue
            values = np.fromstring(line[first_tab + 1 :], sep="\t", dtype=np.float32)
            rows[gene] = values
    if cells is None:
        raise RuntimeError("No header found in count matrix.")
    mat = pd.DataFrame.from_dict(rows, orient="index", columns=cells)
    return mat, cells


def paired_test_table(
    grouped: pd.DataFrame,
    value_cols: list[str],
    group_label: str,
    entity_col: str,
) -> pd.DataFrame:
    out_rows = []
    for entity, sub in grouped.groupby(entity_col, sort=False):
        for col in value_cols:
            wide = sub.pivot_table(
                index="patient_id", columns="treatment_phase", values=col, aggfunc="mean"
            )
            if {"treatment-naive", "post-NACT"}.issubset(wide.columns):
                paired = wide[["treatment-naive", "post-NACT"]].dropna()
            else:
                paired = pd.DataFrame(columns=["treatment-naive", "post-NACT"])
            n_pairs = int(paired.shape[0])
            if n_pairs >= 3:
                diff = paired["post-NACT"] - paired["treatment-naive"]
                mean_naive = paired["treatment-naive"].mean()
                mean_post = paired["post-NACT"].mean()
                mean_diff = diff.mean()
                try:
                    p_t = stats.ttest_1samp(diff, popmean=0.0, nan_policy="omit").pvalue
                except Exception:
                    p_t = np.nan
                try:
                    p_w = stats.wilcoxon(diff).pvalue if np.any(diff != 0) else 1.0
                except Exception:
                    p_w = np.nan
            else:
                mean_naive = np.nan
                mean_post = np.nan
                mean_diff = np.nan
                p_t = np.nan
                p_w = np.nan
            out_rows.append(
                {
                    group_label: entity,
                    "feature": col,
                    "n_pairs": n_pairs,
                    "mean_treatment_naive": mean_naive,
                    "mean_post_NACT": mean_post,
                    "paired_delta_post_minus_naive": mean_diff,
                    "paired_t_p": p_t,
                    "wilcoxon_p": p_w,
                }
            )
    res = pd.DataFrame(out_rows)
    if not res.empty:
        res["paired_t_fdr"] = res.groupby(group_label)["paired_t_p"].transform(bh_adjust)
        res["wilcoxon_fdr"] = res.groupby(group_label)["wilcoxon_p"].transform(bh_adjust)
    return res


def save_heatmap(df: pd.DataFrame, out_file: Path, title: str, center: float = 0.0) -> None:
    plot_df = df.fillna(0)
    vmax = float(np.nanmax(np.abs(plot_df.to_numpy()))) if plot_df.size else 1.0
    vmax = vmax if vmax > 0 else 1.0
    fig, ax = plt.subplots(figsize=(max(7, df.shape[1] * 0.9), max(4, df.shape[0] * 0.35)))
    im = ax.imshow(plot_df.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(plot_df.shape[1]))
    ax.set_xticklabels(plot_df.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(plot_df.shape[0]))
    ax.set_yticklabels(plot_df.index)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("post-NACT minus treatment-naive")
    fig.tight_layout()
    fig.savefig(out_file, dpi=300)
    plt.close(fig)


def scatter_rank_plot(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    hue_col: str,
    size_col: str,
    out_file: Path,
    title: str,
) -> None:
    if df.empty:
        return
    plot_df = df.copy().reset_index(drop=True)
    cats = list(pd.unique(plot_df[hue_col].fillna("NA")))
    cmap = plt.get_cmap("tab20")
    colors = {cat: cmap(i % 20) for i, cat in enumerate(cats)}
    sizes_raw = pd.to_numeric(plot_df[size_col], errors="coerce").fillna(0)
    if sizes_raw.max() > sizes_raw.min():
        sizes = 45 + 180 * (sizes_raw - sizes_raw.min()) / (sizes_raw.max() - sizes_raw.min())
    else:
        sizes = pd.Series(np.full(len(plot_df), 90.0))
    y = np.arange(len(plot_df))
    fig, ax = plt.subplots(figsize=(10, max(5, len(plot_df) * 0.23)))
    for cat in cats:
        mask = plot_df[hue_col].fillna("NA") == cat
        ax.scatter(
            plot_df.loc[mask, x_col],
            y[mask],
            s=sizes.loc[mask],
            color=colors[cat],
            alpha=0.82,
            label=cat,
            edgecolor="white",
            linewidth=0.5,
        )
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df[y_col])
    ax.invert_yaxis()
    ax.set_xlabel(x_col)
    ax.set_title(title)
    ax.legend(loc="best", fontsize=7, frameon=True)
    fig.tight_layout()
    fig.savefig(out_file, dpi=300)
    plt.close(fig)


def df_to_tsv_block(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return "No rows."
    show = df.head(max_rows).copy() if max_rows is not None else df.copy()
    return "```tsv\n" + show.to_csv(sep="\t", index=False) + "```"


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    plt.style.use("ggplot")

    panel = pd.DataFrame(panel_rows())
    pair_df = lr_pairs()
    needed_genes = sorted(set(panel["gene"]) | set(pair_df["ligand"]) | set(pair_df["receptor"]))

    cell_info = pd.read_csv(CELL_INFO, sep="\t")
    cell_info = cell_info.set_index("cell", drop=False)

    count_mat, count_cells = read_selected_counts(COUNTS, needed_genes)
    common_cells = [cell for cell in count_cells if cell in cell_info.index]
    count_mat = count_mat.loc[:, common_cells]
    meta = cell_info.loc[common_cells].copy()

    ncount = meta["nCount_RNA"].to_numpy(dtype=np.float32)
    norm = np.log1p((count_mat.to_numpy(dtype=np.float32) / ncount.reshape(1, -1)) * 1e4)
    norm_df = pd.DataFrame(norm.T, index=common_cells, columns=count_mat.index)
    detected_df = pd.DataFrame((count_mat.to_numpy().T > 0).astype(np.float32), index=common_cells, columns=count_mat.index)

    present = set(count_mat.index)
    panel["present_in_gse165897"] = panel["gene"].isin(present)
    panel.to_csv(TABLES / "protein_coding_panel.tsv", sep="\t", index=False)
    pair_df["ligand_present"] = pair_df["ligand"].isin(present)
    pair_df["receptor_present"] = pair_df["receptor"].isin(present)
    pair_df.to_csv(TABLES / "ligand_receptor_pairs.tsv", sep="\t", index=False)

    audit_scprotrans().to_csv(TABLES / "scprotrans_dependency_audit.tsv", sep="\t", index=False)

    dataset_summary = []
    dataset_summary.append({"metric": "cells", "value": len(meta)})
    dataset_summary.append({"metric": "patients", "value": meta["patient_id"].nunique()})
    dataset_summary.append({"metric": "samples", "value": meta["sample"].nunique()})
    dataset_summary.append({"metric": "panel_genes_requested", "value": len(needed_genes)})
    dataset_summary.append({"metric": "panel_genes_present", "value": len(present)})
    for col in ["treatment_phase", "cell_type", "anatomical_location"]:
        for key, val in meta[col].value_counts().sort_index().items():
            dataset_summary.append({"metric": f"{col}:{key}", "value": int(val)})
    pd.DataFrame(dataset_summary).to_csv(TABLES / "gse165897_dataset_summary.tsv", sep="\t", index=False)

    count_plot = (
        meta.groupby(["treatment_phase", "cell_type"]).size().reset_index(name="n_cells")
    )
    count_plot.to_csv(TABLES / "cell_counts_by_phase_celltype.tsv", sep="\t", index=False)
    count_wide = count_plot.pivot(index="cell_type", columns="treatment_phase", values="n_cells").fillna(0)
    ax = count_wide.plot(kind="bar", figsize=(7.5, 5), color=["#4C78A8", "#F58518"])
    ax.set_title("GSE165897 cells by phase and compartment")
    ax.set_ylabel("n_cells")
    ax.set_xlabel("")
    plt.tight_layout()
    plt.savefig(FIGURES / "fig1_cell_counts_by_phase_celltype.png", dpi=300)
    plt.close()

    group_keys = ["patient_id", "sample", "treatment_phase", "cell_type"]
    gene_group = pd.concat([meta[group_keys], norm_df], axis=1).groupby(group_keys, observed=True).mean().reset_index()
    pct_group = pd.concat([meta[group_keys], detected_df], axis=1).groupby(group_keys, observed=True).mean().reset_index()
    gene_group.to_csv(TABLES / "sample_celltype_gene_logcp10k.tsv.gz", sep="\t", index=False, compression="gzip")
    pct_group.to_csv(TABLES / "sample_celltype_gene_detection_fraction.tsv.gz", sep="\t", index=False, compression="gzip")

    gene_tests = []
    for cell_type, sub in gene_group.groupby("cell_type", observed=True):
        tst = paired_test_table(sub, list(count_mat.index), "cell_type", "cell_type")
        tst = tst[tst["cell_type"] == cell_type].copy()
        gene_tests.append(tst)
    gene_tests_df = pd.concat(gene_tests, ignore_index=True)
    gene_to_module = dict(zip(panel["gene"], panel["module"]))
    gene_tests_df["module"] = gene_tests_df["feature"].map(gene_to_module)
    gene_tests_df.to_csv(TABLES / "paired_gene_delta_by_celltype.tsv", sep="\t", index=False)

    module_scores = meta[group_keys].copy()
    present_panel = panel[panel["present_in_gse165897"]].copy()
    for module, sub in present_panel.groupby("module", sort=False):
        genes = [g for g in sub["gene"] if g in norm_df.columns]
        if genes:
            module_scores[module] = norm_df[genes].mean(axis=1)
    module_group = module_scores.groupby(group_keys, observed=True).mean().reset_index()
    module_group.to_csv(TABLES / "sample_celltype_module_scores.tsv", sep="\t", index=False)

    module_cols = [c for c in module_group.columns if c not in group_keys]
    module_tests = []
    for cell_type, sub in module_group.groupby("cell_type", observed=True):
        tst = paired_test_table(sub, module_cols, "cell_type", "cell_type")
        tst = tst[tst["cell_type"] == cell_type].copy()
        module_tests.append(tst)
    module_tests_df = pd.concat(module_tests, ignore_index=True)
    module_tests_df.to_csv(TABLES / "paired_module_delta_by_celltype.tsv", sep="\t", index=False)

    heat = module_tests_df.pivot(index="feature", columns="cell_type", values="paired_delta_post_minus_naive")
    heat = heat.loc[heat.abs().max(axis=1).sort_values(ascending=False).index]
    save_heatmap(heat, FIGURES / "fig2_module_delta_heatmap.png", "Patient-paired module shifts")

    top_gene = gene_tests_df.dropna(subset=["paired_delta_post_minus_naive"]).copy()
    top_gene["abs_delta"] = top_gene["paired_delta_post_minus_naive"].abs()
    top_gene = top_gene.sort_values(["abs_delta", "n_pairs"], ascending=[False, False]).head(40)
    top_gene.to_csv(TABLES / "top40_gene_delta_candidates.tsv", sep="\t", index=False)
    scatter_rank_plot(
        top_gene,
        "paired_delta_post_minus_naive",
        "feature",
        "cell_type",
        "n_pairs",
        FIGURES / "fig3_top_gene_delta_candidates.png",
        "Top patient-paired gene-level shifts",
    )

    # Ligand-receptor scores use patient-level sender/receiver products.
    lr_rows = []
    cell_types = sorted(gene_group["cell_type"].unique())
    idx_cols = ["patient_id", "treatment_phase", "cell_type"]
    gene_mean = gene_group[idx_cols + list(count_mat.index)].copy()
    for (patient, phase), sub in gene_mean.groupby(["patient_id", "treatment_phase"], observed=True):
        by_type = sub.set_index("cell_type")
        for _, pair in pair_df.iterrows():
            if not (pair["ligand_present"] and pair["receptor_present"]):
                continue
            for sender in cell_types:
                if sender not in by_type.index:
                    continue
                ligand_val = by_type.loc[sender, pair["ligand"]]
                for receiver in cell_types:
                    if receiver not in by_type.index:
                        continue
                    receptor_val = by_type.loc[receiver, pair["receptor"]]
                    lr_rows.append(
                        {
                            "patient_id": patient,
                            "treatment_phase": phase,
                            "sender": sender,
                            "receiver": receiver,
                            "axis": pair["axis"],
                            "ligand": pair["ligand"],
                            "receptor": pair["receptor"],
                            "lr_score": float(ligand_val * receptor_val),
                        }
                    )
    lr_scores = pd.DataFrame(lr_rows)
    lr_scores.to_csv(TABLES / "patient_phase_lr_scores.tsv", sep="\t", index=False)
    lr_scores["interaction"] = (
        lr_scores["sender"]
        + "->"
        + lr_scores["receiver"]
        + ":"
        + lr_scores["ligand"]
        + "-"
        + lr_scores["receptor"]
    )
    lr_tests = paired_test_table(lr_scores, ["lr_score"], "interaction", "interaction")
    if not lr_tests.empty:
        lr_tests["paired_t_fdr_global"] = bh_adjust(lr_tests["paired_t_p"])
        lr_tests["wilcoxon_fdr_global"] = bh_adjust(lr_tests["wilcoxon_p"])
        extra = lr_scores[["interaction", "sender", "receiver", "axis", "ligand", "receptor"]].drop_duplicates()
        lr_tests = lr_tests.merge(extra, on="interaction", how="left")
        lr_tests = lr_tests.sort_values("paired_delta_post_minus_naive", ascending=False)
    lr_tests.to_csv(TABLES / "paired_lr_delta.tsv", sep="\t", index=False)

    top_lr = lr_tests.dropna(subset=["paired_delta_post_minus_naive"]).head(30)
    fig, ax = plt.subplots(figsize=(10, max(5, len(top_lr) * 0.24)))
    y = np.arange(len(top_lr))
    ax.barh(y, top_lr["paired_delta_post_minus_naive"], color="#54A24B")
    ax.set_yticks(y)
    ax.set_yticklabels(top_lr["interaction"])
    ax.invert_yaxis()
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("paired_delta_post_minus_naive")
    ax.set_title("Top post-NACT gains in protein-coding LR proxy scores")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_top_lr_delta.png", dpi=300)
    plt.close(fig)

    # Summarize best candidates by combining gene shifts and LR participation.
    lr_gene_hits = pd.concat(
        [
            lr_tests[["ligand", "paired_delta_post_minus_naive"]].rename(
                columns={"ligand": "feature", "paired_delta_post_minus_naive": "best_lr_delta"}
            ),
            lr_tests[["receptor", "paired_delta_post_minus_naive"]].rename(
                columns={"receptor": "feature", "paired_delta_post_minus_naive": "best_lr_delta"}
            ),
        ],
        ignore_index=True,
    )
    lr_gene_hits = lr_gene_hits.groupby("feature", as_index=False)["best_lr_delta"].max()
    candidates = gene_tests_df.merge(lr_gene_hits, on="feature", how="left")
    candidates["best_lr_delta"] = candidates["best_lr_delta"].fillna(0)
    candidates["direction_bonus"] = np.where(candidates["paired_delta_post_minus_naive"] > 0, 1.0, 0.25)
    candidates["candidate_score"] = (
        candidates["paired_delta_post_minus_naive"].abs().fillna(0)
        * candidates["direction_bonus"]
        * np.log1p(candidates["n_pairs"].fillna(0))
        + np.maximum(candidates["best_lr_delta"], 0)
    )
    candidates = candidates.sort_values("candidate_score", ascending=False)
    candidates.to_csv(TABLES / "ranked_protein_layer_candidates.tsv", sep="\t", index=False)

    top_candidate_plot = candidates.head(35).copy()
    top_candidate_plot["label"] = top_candidate_plot["cell_type"] + ":" + top_candidate_plot["feature"]
    scatter_rank_plot(
        top_candidate_plot,
        "paired_delta_post_minus_naive",
        "label",
        "module",
        "candidate_score",
        FIGURES / "fig5_ranked_candidate_genes.png",
        "Ranked protein-layer candidate genes",
    )

    method_notes = OUT / "scprotrans_hgsoc_independent_summary.md"
    top_modules = module_tests_df.sort_values("paired_delta_post_minus_naive", ascending=False).head(10)
    top_lr_md = lr_tests.head(10) if not lr_tests.empty else pd.DataFrame()
    top_cand_md = candidates.head(15)
    method_notes.write_text(
        "\n".join(
            [
                "# Independent HGSOC scProTrans-ready analysis summary",
                "",
                "## Scope",
                "This run uses GSE165897 public HGSOC scRNA-seq as discovery data.",
                "The official scProTrans repository was cloned and audited, but the cloned",
                "code requires paired RNA/protein data and does not include a pretrained",
                "direct-inference model for HGSOC scRNA alone. Therefore, the tables and",
                "figures produced here are a protein-coding panel / ligand-receptor proxy",
                "baseline designed to be replaced by formal scProTrans predictions once a",
                "trained model or paired reference transfer workflow is available.",
                "",
                "## Dataset",
                f"- Cells: {len(meta)}",
                f"- Patients: {meta['patient_id'].nunique()}",
                f"- Samples: {meta['sample'].nunique()}",
                f"- Present panel genes: {len(present)} / {len(needed_genes)}",
                "",
                "## Top paired module gains",
                df_to_tsv_block(top_modules),
                "",
                "## Top paired ligand-receptor proxy gains",
                df_to_tsv_block(top_lr_md) if not top_lr_md.empty else "No LR tests.",
                "",
                "## Top ranked protein-layer candidates",
                df_to_tsv_block(
                    top_cand_md[
                        [
                            "cell_type",
                            "feature",
                            "module",
                            "n_pairs",
                            "paired_delta_post_minus_naive",
                            "paired_t_p",
                            "paired_t_fdr",
                            "best_lr_delta",
                            "candidate_score",
                        ]
                    ]
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
