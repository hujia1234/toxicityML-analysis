# 网络毒理学联合机器学习与 SHAP 分析工作流

本仓库提供论文研究所用分析流程的源代码，覆盖从基因表达数据预处理、差异表达与加权基因共表达网络分析，到网络毒理学靶点筛选、机器学习建模、SHAP 可解释性分析、免疫浸润及功能富集等环节，便于复现研究结果并理解整体分析思路。

## 分析流程

| 步骤 | 分析内容 | 脚本 |
|:---:|---|---|
| 1 | 探针注释与基因表达矩阵构建 | `04.ann/geo04.ann.pl`、`04.ann/auto_analysis.py` |
| 2 | 表达谱标准化、批次校正与质控 | `geo05.normalize.R`、`geo06.sva.R`、`geo07.boxplot.R`、`geo08.PCA.R` |
| 3 | 差异表达分析 | `geo09.diff.R` |
| 4 | 加权基因共表达网络分析（WGCNA） | `WGCNA_rerun.R` |
| 5 | 疾病相关靶点收集 | `toxicityML11.Disease.R` |
| 6 | 化合物靶点收集与汇总 | `toxicityML16.Compound.R` |
| 7 | 疾病-化合物靶点交集与 Venn 图 | `17.venn/toxicityML17.venn.R`、`17.venn/venn_75genes.R` |
| 8 | GO / KEGG 功能富集分析 | `toxicityML19.GO.R`、`toxicityML20.KEGG.R` |
| 9 | 机器学习特征矩阵构建 | `toxicityML21.MLdata.R` |
| 10 | 多算法机器学习建模与特征频次统计 | `22.ML/refer.ML.R`、`22.ML/toxicityML22.model.R`、`22.ML/gene_frequency.R` |
| 11 | SHAP 可解释性分析 | `23.SHAP/toxicityML24.SHAP.R`、`23.SHAP/plot_TIMP1_dependence.R` |
| 12 | 模型评估与表达核对 | `25.ROC/toxicityML25.ROC.R`、`25.ROC/MMP_diff_check.R` |
| 13 | 免疫细胞浸润分析 | `28.CIBERSORT/SHAP21.CIBERSORT.R`、`28.CIBERSORT/SHAP21.run.R` |
| 14 | GSEA / GSVA 通路分析 | `31GSEA/SHAP19.GSEA.R`、`31GSEA/SHAP21.GSEA.R`、`SHAP20.GSVA.R` |
| 15 | 结果可视化 | `toxicityML24.vol.R`、`SHAP22.barplot.R`、`SHAP23.immuneCor.R` |
| 16 | PPI 网络节点筛选 | `nodeFilter1.R` |

## 目录结构

```
.
├── 04.ann/               探针注释与表达矩阵构建
├── 17.venn/              靶点交集与 Venn 图
├── 22.ML/                多算法机器学习建模
├── 23.SHAP/              SHAP 可解释性分析
├── 25.ROC/               ROC 评估与表达核对
├── 28.CIBERSORT/         免疫细胞浸润分析
├── 31GSEA/               GSEA 富集分析
├── misc/                 辅助脚本
└── *.R                   流程各步骤主脚本
```

## 运行环境

- R 4.4.2
- Python 3.10+
- Perl 5.30+

主要 R 包：`limma`、`sva`、`WGCNA`、`clusterProfiler`、`org.Hs.eg.db`、`enrichplot`、`GSEABase`、`GSVA`、`caret`、`randomForest`、`randomForestSRC`、`glmnet`、`xgboost`、`lightgbm`、`catboost`、`e1071`、`mboost`、`BART`、`Boruta`、`mRMRe`、`ncvreg`、`kernlab`、`pROC`、`DALEX`、`shapviz`、`kernelshap`、`ggplot2`、`ggpubr`、`ggrepel`、`ggvenn`、`ComplexHeatmap`、`pheatmap`、`corrplot`、`dplyr`、`reshape2`。

## 数据与运行说明

- 本仓库仅包含分析代码，不包含表达矩阵、注释表、图表及模型文件等数据与结果文件。
- 各脚本中的输入/输出路径基于研究环境编写，运行前请按实际目录调整。
- 建议按上表顺序执行；同一步骤内的多个脚本按其逻辑先后运行。