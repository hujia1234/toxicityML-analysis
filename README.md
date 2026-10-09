# toxicityML-analysis

网络毒理学联合机器学习与 SHAP 分析的可复现流程代码。

本仓库由 `222/` 与 `555/` 两个研究目录中的分析脚本整理而成，**已移除全部代码注释**，仅保留源码。表达矩阵、注释表、图表、PDF、模型文件等数据与结果文件均未包含在内。

## 目录结构

### 仓库根目录（仅含单个脚本的步骤）

| 脚本 | 说明 |
| --- | --- |
| `geo05.normalize.R` | 表达谱标准化 |
| `geo06.sva.R` | SVA 批次校正 |
| `geo07.boxplot.R` | 标准化箱线图 |
| `geo08.PCA.R` | PCA 分析 |
| `geo09.diff.R` | 差异表达分析 |
| `WGCNA_rerun.R` | WGCNA 加权基因共表达网络分析（重跑版） |
| `toxicityML11.Disease.R` | 疾病靶点整理 |
| `toxicityML16.Compound.R` | 化合物靶点合并与 Venn |
| `toxicityML19.GO.R` | GO 富集分析 |
| `toxicityML20.KEGG.R` | KEGG 通路富集分析 |
| `toxicityML21.MLdata.R` | 机器学习数据准备 |
| `toxicityML24.vol.R` | 火山图 |
| `SHAP22.barplot.R` | 免疫浸润柱状图 |
| `SHAP23.immuneCor.R` | 核心基因与免疫细胞相关性 |
| `SHAP20.GSVA.R` | GSVA 富集分析 |
| `nodeFilter1.R` | PPI 网络节点筛选 |

### 子目录（含多个脚本的步骤）

| 目录 | 脚本 | 说明 |
| --- | --- | --- |
| `04.ann/` | `geo04.ann.pl`、`auto_analysis.py` | 探针注释与基因表达矩阵构建 |
| `17.venn/` | `toxicityML17.venn.R`、`venn_75genes.R` | 交集基因 Venn 图 |
| `22.ML/` | `refer.ML.R`、`toxicityML22.model.R`、`gene_frequency.R` | 多算法机器学习建模（优化版） |
| `23.SHAP/` | `toxicityML24.SHAP.R`、`plot_TIMP1_dependence.R` | SHAP 可解释性分析（共识版） |
| `25.ROC/` | `toxicityML25.ROC.R`、`MMP_diff_check.R` | ROC 曲线与基因表达核对 |
| `28.CIBERSORT/` | `SHAP21.CIBERSORT.R`、`SHAP21.run.R` | 免疫细胞浸润分析 |
| `31GSEA/` | `SHAP19.GSEA.R`、`SHAP21.GSEA.R` | GSEA 富集分析 |
| `misc/` | `_check_module_targets*.R`、`_mmgs.R`、`_extract_pdfs.py` | 辅助与调试脚本 |

## 整理说明

- **两套目录合并**：`222/` 为完整流程（编号 01–32），`555/` 为其中若干步骤的重跑/优化版本。合并后同一编号只保留一份，优先采用 `555/` 的重跑/优化版本（如 `WGCNA_rerun.R`、`22.ML`、`23.SHAP`、`SHAP20.GSVA.R`）。
- **取消单文件目录**：仅含一个脚本的步骤目录已取消，脚本直接放在仓库根目录；含多个脚本的步骤目录保留。
- **已排除**：
  - `AutoDockGUI/`：独立的分子对接桌面工具（Python GUI），不属于本分析流程。
  - `toxicityML13.ChEMBL.pl`：PAR 打包的 Perl 可执行文件（二进制，非源码）。
  - `04.ann/GSE38713 - 副本/geo04.ann.pl`：备份副本，且内嵌了大段注释数据。
  - 各类数据/结果文件（`.txt`、`.pdf`、`.png`、`.xls`、`.rds` 等）。
- **注释处理**：删除 R / Perl / Python 中的 `#` 注释、Perl 的 POD 文档块与 Python 文档字符串；保留 shebang 行。所有脚本在删除注释后均已通过语法校验（R `parse`、Python `compile`、Perl `-c`）。

## 环境

- R 4.4.2
- Python 3.10+
- Perl 5.30+