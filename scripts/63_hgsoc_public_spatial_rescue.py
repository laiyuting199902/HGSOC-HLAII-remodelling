#!/usr/bin/env python3

"""公共 HGSOC 空间队列的 HLA-II 可行性与边界审计。"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread


CORE_GENES = ["CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1"]
EPITHELIAL_GENES = ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "PAX8", "WFDC2", "MUC16", "MSLN"]
IMMUNE_GENES = ["PTPRC", "CD3D", "CD3E", "CD8A", "CD8B", "LST1", "TYROBP", "AIF1", "MS4A1", "CD79A"]
MODULES = {
    "epithelial": EPITHELIAL_GENES,
    "t_cd8": ["CD3D", "CD3E", "CD8A", "CD8B", "CCL5", "NKG7"],
    "b_tls": ["MS4A1", "CD79A", "CD79B", "CD37", "CXCL13", "CCL19", "LTB"],
    "myeloid": ["AIF1", "LST1", "TYROBP", "FCER1G", "CTSS", "C1QA", "C1QB", "C1QC"],
}


def parse_args() -> argparse.Namespace:
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("data/raw/spatial_public"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=project / "outputs" / "scprotrans_hgsoc_v4" / "spatial_public_rescue",
    )
    return parser.parse_args()


def log_normalize(matrix: sparse.spmatrix) -> sparse.csc_matrix:
    matrix = matrix.tocsc().astype(np.float64)
    library = np.asarray(matrix.sum(axis=0)).ravel()
    scale = 10000.0 / np.maximum(library, 1.0)
    normalized = matrix @ sparse.diags(scale)
    normalized.data = np.log1p(normalized.data)
    return normalized


def module_score(matrix: sparse.csc_matrix, genes: list[str], members: list[str]) -> np.ndarray:
    lookup = {gene: index for index, gene in enumerate(genes)}
    indexes = [lookup[gene] for gene in members if gene in lookup]
    if not indexes:
        return np.full(matrix.shape[1], np.nan)
    return np.asarray(matrix[indexes, :].mean(axis=0)).ravel()


def read_features(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        rows = [line.rstrip().split("\t") for line in handle]
    return [row[1] if len(row) > 1 else row[0] for row in rows]


def visium_sample_scores(feature_path: Path, matrix_path: Path) -> dict[str, np.ndarray | int | bool]:
    genes = read_features(feature_path)
    matrix = mmread(matrix_path).tocsc()
    if matrix.shape[0] != len(genes):
        raise RuntimeError(f"特征数与矩阵维度不一致：{matrix_path}")
    normalized = log_normalize(matrix)
    core = module_score(normalized, genes, CORE_GENES)
    epithelial = module_score(normalized, genes, EPITHELIAL_GENES)
    immune = module_score(normalized, genes, IMMUNE_GENES)
    purity_proxy = epithelial - immune
    tumor_proxy = purity_proxy >= np.nanquantile(purity_proxy, 0.70)
    return {
        "n_spots": matrix.shape[1],
        "n_tumor_proxy_spots": int(tumor_proxy.sum()),
        "all_core_mean": float(np.nanmean(core)),
        "tumor_proxy_core_mean": float(np.nanmean(core[tumor_proxy])),
        "tumor_proxy_epithelial_mean": float(np.nanmean(epithelial[tumor_proxy])),
        "tumor_proxy_immune_mean": float(np.nanmean(immune[tumor_proxy])),
        "tumor_proxy_core_detected_fraction": float(np.mean(core[tumor_proxy] > 0)),
        "all_core_genes_covered": all(gene in set(genes) for gene in CORE_GENES),
    }


def analyze_gse189843(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for number in range(485, 497):
        gsm = f"GSM5708{number}"
        feature_path = next(root.glob(f"{gsm}_features_*.tsv.gz"))
        matrix_path = next(root.glob(f"{gsm}_matrix_*.mtx.gz"))
        values = visium_sample_scores(feature_path, matrix_path)
        rows.append({"cohort": "GSE189843", "sample": gsm, "stage": "chemo-naive", "response": "excellent" if number <= 490 else "poor", **values})
    samples = pd.DataFrame(rows)
    excellent = samples.loc[samples["response"].eq("excellent"), "tumor_proxy_core_mean"]
    poor = samples.loc[samples["response"].eq("poor"), "tumor_proxy_core_mean"]
    test = stats.mannwhitneyu(excellent, poor, alternative="two-sided")
    tests = pd.DataFrame([
        {
            "cohort": "GSE189843",
            "comparison": "excellent_vs_poor_NACT_response",
            "unit": "patient/sample",
            "n_group_1": len(excellent),
            "n_group_2": len(poor),
            "mean_group_1": excellent.mean(),
            "mean_group_2": poor.mean(),
            "delta_group_1_minus_group_2": excellent.mean() - poor.mean(),
            "statistic": test.statistic,
            "p_value": test.pvalue,
            "interpretation": "治疗前空间 RNA 疗效分层未通过；保留为阴性结果",
        }
    ])
    return samples, tests


def analyze_gse211956(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    responses = ["poor", "good", "good", "partial", "good", "partial", "poor", "poor"]
    rows = []
    for spatial_id in range(1, 9):
        prefix = f"GSM{6506109 + spatial_id}_SP{spatial_id}"
        values = visium_sample_scores(root / f"{prefix}_features.tsv.gz", root / f"{prefix}_matrix.mtx.gz")
        rows.append({"cohort": "GSE211956", "sample": f"SP{spatial_id}", "stage": "IDS", "response": responses[spatial_id - 1], **values})
    samples = pd.DataFrame(rows)
    poor = samples.loc[samples["response"].eq("poor"), "tumor_proxy_core_mean"]
    nonpoor = samples.loc[~samples["response"].eq("poor"), "tumor_proxy_core_mean"]
    test = stats.mannwhitneyu(poor, nonpoor, alternative="two-sided")
    tests = pd.DataFrame([
        {
            "cohort": "GSE211956",
            "comparison": "poor_vs_good_or_partial_NACT_response",
            "unit": "patient/sample",
            "n_group_1": len(poor),
            "n_group_2": len(nonpoor),
            "mean_group_1": poor.mean(),
            "mean_group_2": nonpoor.mean(),
            "delta_group_1_minus_group_2": poor.mean() - nonpoor.mean(),
            "statistic": test.statistic,
            "p_value": test.pvalue,
            "interpretation": "IDS 空间 RNA 中差疗效组较低，但样本量不足且未通过检验",
        }
    ])
    return samples, tests


def parse_geomx_site(sample_site: str) -> str:
    for site in ["Fallopian Tube", "Omentum", "Ovary", "Lymph Node"]:
        if sample_site.endswith(site):
            return site
    return "Other"


def analyze_gse276935(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    segments = pd.read_excel(path, sheet_name="SegmentProperties")
    matrix = pd.read_excel(path, sheet_name="TargetCountMatrix").set_index("TargetName")
    missing = [gene for gene in CORE_GENES if gene not in matrix.index]
    if missing:
        raise RuntimeError(f"GSE276935 缺少核心基因：{', '.join(missing)}")
    expression = np.log2(matrix + 1.0).T
    score_sets = {"hlaii_core": CORE_GENES, **MODULES}
    for score_name, members in score_sets.items():
        covered = [gene for gene in members if gene in expression.columns]
        expression[score_name] = expression[covered].mean(axis=1)
    metadata = segments.set_index("TargetName").loc[expression.index].copy()
    metadata = metadata.join(expression[list(score_sets)])
    metadata["sample_site"] = metadata.index.to_series().str.replace(
        r"\d{3}(CD20|Immune|Stroma|Tumor|Full ROI|Other)$", "", regex=True
    )
    metadata["site"] = metadata["sample_site"].map(parse_geomx_site)
    metadata["roi_key"] = metadata["sample_site"] + "_" + metadata["ROILabel"].astype(str)
    metadata["qc_clean"] = metadata["QCFlags"].isna()
    output_columns = [
        "sample_site", "site", "ROILabel", "roi_key", "SegmentLabel", "qc_clean",
        "hlaii_core", *CORE_GENES, "epithelial", "t_cd8", "b_tls", "myeloid",
    ]
    segment_scores = metadata.join(expression[CORE_GENES])[output_columns].reset_index(names="aoi")

    clean_metadata = metadata.loc[metadata["qc_clean"]].copy()
    sample_scores = (
        clean_metadata.join(expression[CORE_GENES])
        .groupby(["sample_site", "site", "SegmentLabel"], observed=True)[list(score_sets) + CORE_GENES]
        .mean()
        .reset_index()
    )
    tumor = sample_scores.loc[sample_scores["SegmentLabel"].eq("Tumor")].copy()

    roi_wide = clean_metadata.pivot_table(
        index=["sample_site", "roi_key"], columns="SegmentLabel",
        values=["hlaii_core", "epithelial", "t_cd8", "b_tls", "myeloid"], aggfunc="mean",
    )
    association_rows = []
    for neighbor_segment in ["Immune", "CD20", "Stroma"]:
        for neighbor_module in ["t_cd8", "b_tls", "myeloid"]:
            tumor_values = roi_wide.get(("hlaii_core", "Tumor"))
            neighbor_values = roi_wide.get((neighbor_module, neighbor_segment))
            if tumor_values is None or neighbor_values is None:
                continue
            paired = pd.concat([tumor_values, neighbor_values], axis=1).dropna().reset_index()
            patient_equal = paired.groupby("sample_site")[[paired.columns[-2], paired.columns[-1]]].mean()
            if len(patient_equal) >= 4:
                rho, p_value = stats.spearmanr(patient_equal.iloc[:, 0], patient_equal.iloc[:, 1])
            else:
                rho, p_value = np.nan, np.nan
            association_rows.append(
                {
                    "cohort": "GSE276935",
                    "tumor_feature": "PanCK-segmented tumor HLA-II core",
                    "neighbor_segment": neighbor_segment,
                    "neighbor_module": neighbor_module,
                    "n_samples": len(patient_equal),
                    "spearman_rho": rho,
                    "p_value": p_value,
                    "inference_unit": "sample_site_equal",
                }
            )
    associations = pd.DataFrame(association_rows)
    associations["fdr"] = benjamini_hochberg(associations["p_value"].to_numpy())

    tumor_aoi = clean_metadata.loc[clean_metadata["SegmentLabel"].eq("Tumor")]
    aoi_rho, aoi_p_value = stats.spearmanr(tumor_aoi["hlaii_core"], tumor_aoi["epithelial"])
    tumor_sample = tumor[["sample_site", "hlaii_core", "epithelial"]].dropna()
    sample_rho, sample_p_value = stats.spearmanr(tumor_sample["hlaii_core"], tumor_sample["epithelial"])
    all_tumor_aoi = metadata.loc[metadata["SegmentLabel"].eq("Tumor")]
    all_aoi_rho, all_aoi_p_value = stats.spearmanr(all_tumor_aoi["hlaii_core"], all_tumor_aoi["epithelial"])
    audit = pd.DataFrame([
        {
            "cohort": "GSE276935",
            "analysis": "tumor_AOI_HLAII_vs_epithelial_module",
            "unit": "QC_passing_AOI_exploratory",
            "n": len(tumor_aoi),
            "estimate": aoi_rho,
            "p_value": aoi_p_value,
            "claim_boundary": "QC 合格 PanCK 肿瘤 AOI 的探索性结果；不能作为患者级独立验证或蛋白证据",
        },
        {
            "cohort": "GSE276935",
            "analysis": "tumor_sample_HLAII_vs_epithelial_module",
            "unit": "QC_passing_sample_site_equal",
            "n": len(tumor_sample),
            "estimate": sample_rho,
            "p_value": sample_p_value,
            "claim_boundary": "仅纳入至少一个 QC 合格肿瘤 AOI 的 sample-site，并按 sample-site 等权",
        },
        {
            "cohort": "GSE276935",
            "analysis": "PanCK_tumor_sample_coverage",
            "unit": "QC_passing_sample_site",
            "n": tumor["sample_site"].nunique(),
            "estimate": tumor["hlaii_core"].mean(),
            "p_value": np.nan,
            "claim_boundary": "仅统计至少一个 QC 合格 PanCK 肿瘤 AOI 且五基因完整覆盖的 sample-site",
        },
        {
            "cohort": "GSE276935",
            "analysis": "all_tumor_AOI_sensitivity",
            "unit": "all_AOI_sensitivity",
            "n": len(all_tumor_aoi),
            "estimate": all_aoi_rho,
            "p_value": all_aoi_p_value,
            "claim_boundary": "包含 QC 标记 AOI 的敏感性结果，不用于主要推断",
        },
    ])
    return segment_scores, sample_scores, associations, audit


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    finite = np.isfinite(p_values)
    result = np.full_like(p_values, np.nan)
    if not finite.any():
        return result
    valid = p_values[finite]
    order = np.argsort(valid)
    ranked = valid[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    valid_result = np.empty_like(adjusted)
    valid_result[order] = np.minimum(adjusted, 1.0)
    result[finite] = valid_result
    return result


def write_registry(output_dir: Path) -> None:
    registry = pd.DataFrame(
        [
            ["GSE189843", "Visium RNA", 12, "chemo-naive", False, True, False, "accessible", "肿瘤富集 spot 的 HLA-II RNA 定位与 NACT 疗效背景"],
            ["GSE211956", "Visium RNA", 8, "IDS", False, True, False, "accessible", "IDS 肿瘤富集 spot 的 HLA-II RNA 与空间邻域"],
            ["GSE276935", "GeoMx WTA RNA", 11, "cross-sectional", False, True, False, "accessible", "PanCK 分割肿瘤区的 HLA-II RNA 和同 ROI 免疫区"],
            ["syn66694443", "t-CyCIF single-cell spatial protein", 16, "primary/interval", True, True, True, "downloaded_analyzed", "项目 syn53283672；6 对患者可配对；广义肿瘤 MHC-II 5/6 上升"],
            ["syn72380119", "single-cell spatial protein", 280, "cross-sectional", False, True, True, "public_login_required", "929 张空间图；直接肿瘤 MHC-II 蛋白和临床结局"],
        ],
        columns=[
            "accession", "platform", "reported_patient_n", "treatment_context",
            "paired_pre_post", "tumor_compartment_resolvable", "measured_mhcii_protein",
            "access_status", "recommended_role",
        ],
    )
    registry.to_csv(output_dir / "public_spatial_dataset_registry.tsv", sep="\t", index=False)


def write_chinese_report(
    output_dir: Path,
    visium_samples: pd.DataFrame,
    response_tests: pd.DataFrame,
    geomx_audit: pd.DataFrame,
    associations: pd.DataFrame,
) -> None:
    gse189 = response_tests.loc[response_tests["cohort"].eq("GSE189843")].iloc[0]
    gse211 = response_tests.loc[response_tests["cohort"].eq("GSE211956")].iloc[0]
    aoi = geomx_audit.loc[geomx_audit["analysis"].eq("tumor_AOI_HLAII_vs_epithelial_module")].iloc[0]
    patient = geomx_audit.loc[geomx_audit["analysis"].eq("tumor_sample_HLAII_vs_epithelial_module")].iloc[0]
    finite_associations = associations.loc[associations["p_value"].notna()].sort_values("p_value")
    strongest = finite_associations.iloc[0] if not finite_associations.empty else None
    adjacency_text = (
        f"QC 合格 AOI 的相邻分区审计样本量很小；可检验组合中最小 P 值为 "
        f"{strongest['p_value']:.3f}（FDR = {strongest['fdr']:.3f}，n = {int(strongest['n_samples'])}），"
        "未形成可晋级的空间生态关联。"
        if strongest is not None
        else "QC 合格 AOI 的相邻分区配对不足，无法进行稳定的空间生态关联检验。"
    )
    text = f"""# 公共空间数据分析：可行性与证据边界

## 结论

公共空间数据提供三类相互区分的证据：

1. `GSE276935`、`GSE189843` 和 `GSE211956` 能支持 **HLA-II 的空间转录定位**。
2. `syn53283672` 项目中的更新文件 `syn66694443` 已下载并完成分析，能够提供 **6 对患者的纵向肿瘤细胞 MHC-II 直接空间蛋白测量**；`syn72380119` 仍需单独申请访问。
3. 公共观察性空间数据不能替代抗原呈递功能实验，也不能证明 HLA-II 导致免疫激活。

## 已完成的原始数据审计

- 三个无需登录的公开队列均完整覆盖 CD74、HLA-DRA、HLA-DRB1、HLA-DPA1 和 HLA-DPB1。
- `GSE276935` 含 566 个 GeoMx AOI，其中 41 个为 PanCK 分割的肿瘤 AOI。主要分析严格排除带有厂商 QC 标记的 AOI，最终保留 {int(aoi['n'])} 个 QC 合格肿瘤 AOI、覆盖 {int(patient['n'])} 个 sample-site。HLA-II 核心与独立上皮模块在 AOI 层面呈正相关（rho = {aoi['estimate']:.3f}，P = {aoi['p_value']:.4g}），按 sample-site 等权后仍为正（rho = {patient['estimate']:.3f}，P = {patient['p_value']:.3f}，n = {int(patient['n'])}）。主要价值是 PanCK 肿瘤分区中的转录定位；AOI 层面 P 值仅为探索性，患者等权结果才是主要推断单位。
- 预设审计了 3 类相邻分区和 3 个免疫模块。{adjacency_text}
- `GSE189843` 的 12 例治疗前 Visium 中，NACT 优/差疗效组的肿瘤富集 spot HLA-II 核心无显著差异（P = {gse189['p_value']:.3f}）。
- `GSE211956` 的 8 例 IDS Visium 中，差疗效组平均 HLA-II 较低，但未通过患者级检验（P = {gse211['p_value']:.3f}）。因此空间数据目前不能把五基因核心晋级为疗效预测标志物。
- `syn66694443` 含 6,422,542 个细胞和 16 例患者，其中 6 例同时具有 primary 与 interval 样本。广义肿瘤区室 MHC-II 为 5/6 例上升，患者等权平均差值 0.084（t 分布 95% CI -0.002 至 0.170；t 检验 P = 0.0535；精确 Wilcoxon P = 0.0938）。严格上皮敏感性同样为 5/6 例上升，平均差值 0.097（95% CI 0.010 至 0.185；t 检验 P = 0.0359；精确 Wilcoxon P = 0.0625）。

## 纵向空间蛋白证据

`syn66694443` 与 GSE266577 来自同一 Cancer Cell 研究，但提供不同测量技术：高通量 t-CyCIF 直接测量 MHC-II 蛋白、保留单细胞坐标和肿瘤/免疫/基质身份，并有 6 个可分析的患者配对。分析检验：

- 肿瘤细胞 MHC-II 蛋白是否在 IDS 后升高；
- 结果在同部位配对和所有配对中是否一致；
- MHC-II 高肿瘤细胞是否靠近 CD8 T 细胞、巨噬细胞或肿瘤-基质界面；
- RNA 五基因变化与同患者 t-CyCIF MHC-II 蛋白变化是否一致。

这些数据连接了纵向肿瘤细胞 RNA 与纵向空间实测 MHC-II 蛋白。但 6 对蛋白样本与 RNA 队列来自同一研究且患者重叠，因此属于不同测量技术支持，不是独立队列验证；空间邻近仍不是功能验证。

## 结果解释边界

- 可写：`公共空间转录数据支持 CD74/HLA-II 信号存在于 PanCK 界定的 HGSOC 肿瘤区，并与独立上皮转录背景一致。`
- 可写：`患者配对空间蛋白成像直接评估了肿瘤细胞 MHC-II 蛋白的治疗相关变化，并提供方向一致的支持。`
- 不可写：`空间转录证明 HLA-II 蛋白升高`、`空间共定位证明抗原呈递功能` 或 `HLA-II 导致化疗反应。`

## 数据质量说明

- 主要 GeoMx 结果只纳入 `QCFlags` 为空的 AOI；包含全部 41 个肿瘤 AOI 的结果仅作为敏感性分析。
- Visium 的“肿瘤富集 spot”由与 HLA-II 终点无关的上皮减免疫转录代理评分界定，并不等同于病理学逐 spot 肿瘤标注。
- 所有显著性检验均以患者或 sample-site 为推断单位；AOI 和 spot 不作为独立患者重复。
"""
    (output_dir / "公共空间数据分析评估.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    gse189_samples, gse189_tests = analyze_gse189843(args.data_root / "GSE189843")
    gse211_samples, gse211_tests = analyze_gse211956(args.data_root / "GSE211956")
    visium_samples = pd.concat([gse189_samples, gse211_samples], ignore_index=True)
    response_tests = pd.concat([gse189_tests, gse211_tests], ignore_index=True)
    visium_samples.to_csv(args.output_dir / "visium_hlaii_sample_summary.tsv", sep="\t", index=False)
    response_tests.to_csv(args.output_dir / "visium_response_tests.tsv", sep="\t", index=False)

    geomx_segments, geomx_samples, associations, geomx_audit = analyze_gse276935(
        args.data_root / "GSE276935_Processed_Data_File_All_Experiments.xlsx"
    )
    geomx_segments.to_csv(args.output_dir / "gse276935_geomx_aoi_scores.tsv", sep="\t", index=False)
    geomx_samples.to_csv(args.output_dir / "gse276935_geomx_sample_segment_scores.tsv", sep="\t", index=False)
    associations.to_csv(args.output_dir / "gse276935_patient_equal_spatial_associations.tsv", sep="\t", index=False)
    geomx_audit.to_csv(args.output_dir / "gse276935_hlaii_audit.tsv", sep="\t", index=False)
    write_registry(args.output_dir)
    write_chinese_report(args.output_dir, visium_samples, response_tests, geomx_audit, associations)
    print(f"空间数据分析完成：{args.output_dir}")


if __name__ == "__main__":
    main()
