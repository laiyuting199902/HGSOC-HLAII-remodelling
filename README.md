# HGSOC 新辅助化疗后 CD74/HLA-II 重塑分析代码

本仓库仅保存公共数据分析所需的 R、Python 和 C++ 代码，以及运行环境、测试和数据获取说明。原始数据、分析结果、图件、稿件及投稿文件均不在本仓库中分发。

## 目录

- `scripts/`：单细胞、整体组织、空间转录组、蛋白组和 scProTrans 迁移分析脚本。
- `R/`：可复用的 R 辅助函数。
- `tools/`：大规模稀疏矩阵流式处理工具。
- `tests/`：无需原始数据的单元测试，以及需要完整公共数据的集成测试。
- `workflow/`：主要分析、外部数据分析和测试入口。
- `environment/`：Python、R 及相关依赖版本。
- `manifests/public_data_accessions.tsv`：公共数据登录号和获取地址。
- `manifests/analysis_workflow.tsv`：分析步骤与脚本顺序。

## 数据获取

本项目未产生新的患者数据。请根据 `manifests/public_data_accessions.tsv` 从 GEO、Synapse、GDC、cBioPortal/CPTAC 或 Human Protein Atlas 下载原始数据，并放入 `data/raw/`。受 Synapse 条款管理的文件应由每位使用者登录后自行获取。

建议目录结构：

```text
data/raw/
  gse266577/
  gse165897/
  gse201047/
  gse227666/
  gse319500/
  gse143897/
  spatial_public/
  spatial_protein/syn66694443/
  scprotrans_reference/
  msigdb/2026.1.Hs/
  tcga_ov/
  cbioportal_cptac/
```

## 运行

主要配对单细胞分析：

```bash
bash workflow/run_primary.sh
```

外部队列、空间和实测蛋白分析：

```bash
bash workflow/run_external.sh
```

查看 scProTrans 迁移分析的命令接口：

```bash
bash workflow/run_scprotrans_interfaces.sh
```

脚本生成的 `outputs/`、`results/`、`reports/` 和 `figures/` 均由 Git 忽略。

## 测试

```bash
bash workflow/run_unit_tests.sh
```

`tests/integration/` 中的检查需要先运行相应分析并生成本地结果。

## 第三方实现

scProTrans/ProTrans 原始实现来自 <https://github.com/MengyuanZhaoo/ProTrans>。本分析使用的参考提交为 `65bfb3c52a7261f5e5c7e64246c6d1b276805e2a`，第三方源码未复制进本仓库。

## 许可证

本仓库自编分析代码按 MIT License 发布。公共数据仍受各来源数据库及原始研究的使用条款约束。
