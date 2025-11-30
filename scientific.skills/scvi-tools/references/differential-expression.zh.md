<!-- 此文件由机器翻译自 differential-expression.md -->

# scvi-tools 中的差异表达分析

本文档提供了有关使用 scvi-tools 概率框架进行差异表达 (DE) 分析的详细信息。

## 概述

scvi-tools 实现贝叶斯差异表达测试，利用学习的生成模型来估计组之间的表达差异。与传统方法相比，这种方法具有以下几个优点：

- **批量校正**：对批量校正表示进行 DE 测试
- **不确定性量化**：效应大小的概率估计
- **零通货膨胀处理**：辍学和零的正确建模
- **灵活的比较**：任何组或细胞类型之间
- **多种模式**：适用于 RNA、蛋白质 (totalVI) 和可及性 (PeakVI)

## 核心统计框架

### 问题定义

目标是估计两个条件之间表达的对数倍数变化：

```
log fold-change = log(μ_B) - log(μ_A)
```

其中 μ_A 和 μ_B 是条件 A 和 B 下的平均表达水平。

### 三阶段过程

**第 1 阶段：估计表达水平**
- 细胞状态后验分布的样本
- 从学习的生成模型生成表达值
- 跨细胞聚合以获得群体水平的估计

**第 2 阶段：检测相关特征（假设检验）**
- 使用贝叶斯框架测试差异表达
- 两种测试模式可供选择：
  - **“普通”模式**：点原假设 (β = 0)
  - **“改变”模式**：复合假设 (|β| ≤ δ)

**第三阶段：控制错误发现**
- 后验预期错误发现比例（FDP）控制
- 选择最大发现数，确保 E[FDP] ≤ α

## 基本用法

### 简单的两组比较

<<<代码块_1>>>

### 一与其他比较

<<<代码块_2>>>

### 所有成对比较

<<<代码块_3>>>

## 关键参数

### `groupby`（必需）
`adata.obs` 中的列定义要比较的组。

<<<代码块_4>>>

### `group1` 和 `group2`
进行分组比较。如果 `group2` 为 None，则将 `group1` 与所有其他值进行比较。

<<<代码块_5>>>

### `mode`（假设检验模式）

**“vanilla”模式**（默认）：点原假设
- 测试 β 是否准确地 = 0
- 更敏感，但可能会发现微不足道的小影响

**“改变”模式**：复合原假设
- 测试是否|β| ≤δ
- 需要具有生物学意义的改变
- 减少微小效应的错误发现

<<<代码块_6>>>

### `delta`
“更改”模式的最小效果大小阈值。
- 典型值：0.25、0.5、0.7（对数刻度）
- log2(1.5) ≈ 0.58（1.5 倍变化）
- log2(2) = 1.0（2 倍变化）

```python
# Require at least 1.5-fold change
de = model.differential_expression(
    groupby="condition",
    group1="disease",
    group2="healthy",
    mode="change",
    delta=0.58  # log2(1.5)
)
```

### `fdr_target`
错误发现率阈值（默认：0.05）

```python
# More stringent FDR control
de = model.differential_expression(
    groupby="cell_type",
    group1="T cells",
    fdr_target=0.01
)
```

### `batch_correction`
DE测试时是否进行批量校正（默认：True）

```python
# Test within a specific batch
de = model.differential_expression(
    groupby="cell_type",
    group1="T cells",
    group2="B cells",
    batch_correction=False
)
```

### `n_samples`
用于估计的后验样本数量（默认值：5000）
- 更多样本=更准确但更慢
- 降低速度，提高精度

```python
# High precision analysis
de = model.differential_expression(
    groupby="cell_type",
    group1="T cells",
    n_samples=10000
)
```

## 解释结果

### 输出列

结果 DataFrame 包含几个重要的列：

**效果大小估计**：
- `lfc_mean`：平均对数倍数变化
- `lfc_median`：中值对数倍数变化
- `lfc_std`：对数倍数变化的标准偏差
- `lfc_min`：效果大小的下限
- `lfc_max`：效果大小的上限

**统计意义**：
- `bayes_factor`：差异表达的贝叶斯因子
  - 更高的值=更有力的证据
  - >3 通常被认为有意义
- `is_de_fdr_0.05`：布尔值，指示基因是否在 FDR 0.05 处为 DE
- `is_de_fdr_0.1`：布尔值，指示基因在 FDR 0.1 处是否为 DE

**表达水平**：
- `mean1`：组 1 中的平均表达式
- `mean2`：组 2 中的平均表达式
- `non_zeros_proportion1`：组 1 中非零单元格的比例
- `non_zeros_proportion2`：第 2 组中非零单元格的比例

### 示例解释

```python
de_results = model.differential_expression(
    groupby="cell_type",
    group1="T cells",
    group2="B cells"
)

# Find significantly upregulated genes in T cells
upreg_tcells = de_results[
    (de_results["is_de_fdr_0.05"]) &
    (de_results["lfc_mean"] > 0)
].sort_values("lfc_mean", ascending=False)

print(f"Upregulated genes in T cells: {len(upreg_tcells)}")
print(upreg_tcells.head(10))

# Find genes with large effect sizes
large_effect = de_results[
    (de_results["is_de_fdr_0.05"]) &
    (abs(de_results["lfc_mean"]) > 1)  # 2-fold change
]
```

## 高级用法

### 特定单元格内的 DE

```python
# Test DE only within a subset of cells
subset_indices = adata.obs["tissue"] == "lung"

de = model.differential_expression(
    idx1=adata.obs["cell_type"] == "T cells" & subset_indices,
    idx2=adata.obs["cell_type"] == "B cells" & subset_indices
)
```

### 特定批次的 DE

```python
# Test DE within each batch separately
batches = adata.obs["batch"].unique()

batch_de_results = {}
for batch in batches:
    batch_idx = adata.obs["batch"] == batch
    batch_de_results[batch] = model.differential_expression(
        idx1=(adata.obs["condition"] == "treated") & batch_idx,
        idx2=(adata.obs["condition"] == "control") & batch_idx
    )
```

### 伪散装 DE

```python
# Aggregate cells before DE testing
# Useful for low cell counts per group

de = model.differential_expression(
    groupby="cell_type",
    group1="rare_cell_type",
    group2="common_cell_type",
    n_samples=10000,  # More samples for stability
    batch_correction=True
)
```

## 可视化

### 火山图

```python
import matplotlib.pyplot as plt
import numpy as np

de = model.differential_expression(
    groupby="condition",
    group1="treated",
    group2="control"
)

# Volcano plot
plt.figure(figsize=(10, 6))
plt.scatter(
    de["lfc_mean"],
    -np.log10(1 / (de["bayes_factor"] + 1)),
    c=de["is_de_fdr_0.05"],
    cmap="coolwarm",
    alpha=0.5
)
plt.xlabel("Log Fold Change")
plt.ylabel("-log10(1/Bayes Factor)")
plt.title("Volcano Plot: Treated vs Control")
plt.axvline(x=0, color='k', linestyle='--', linewidth=0.5)
plt.show()
```

### 顶级 DE 基因热图

```python
import seaborn as sns

# Get top DE genes
top_genes = de.sort_values("lfc_mean", ascending=False).head(50).index

# Get normalized expression
norm_expr = model.get_normalized_expression(
    adata,
    indices=adata.obs["condition"].isin(["treated", "control"]),
    gene_list=top_genes
)

# Plot heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(
    norm_expr.T,
    cmap="viridis",
    xticklabels=False,
    yticklabels=top_genes
)
plt.title("Top 50 DE Genes")
plt.show()
```

### 排名基因图

```python
# Plot genes ranked by effect size
de_sorted = de.sort_values("lfc_mean", ascending=False)

plt.figure(figsize=(12, 6))
plt.plot(range(len(de_sorted)), de_sorted["lfc_mean"].values)
plt.axhline(y=0, color='r', linestyle='--')
plt.xlabel("Gene Rank")
plt.ylabel("Log Fold Change")
plt.title("Genes Ranked by Effect Size")
plt.show()
```

## 与传统方法的比较

### scvi-tools 与 Wilcoxon 测试
```python
import scanpy as sc

# Traditional Wilcoxon test
sc.tl.rank_genes_groups(
    adata,
    groupby="cell_type",
    method="wilcoxon",
    key_added="wilcoxon"
)

# scvi-tools DE
de_scvi = model.differential_expression(
    groupby="cell_type",
    group1="T cells"
)

# Compare results
wilcox_results = sc.get.rank_genes_groups_df(adata, group="T cells", key="wilcoxon")
```

**scvi-tools 的优点**：
- 自动考虑批次效果
- 正确处理零通胀
- 提供不确定性量化
- 不需要任意伪计数
- 更好的统计特性

**何时使用 Wilcoxon**：
- 非常快速的探索性分析
- 使用 Wilcoxon 与已发表的结果进行比较

## 多模式 DE

### 蛋白质 DE（总 VI）

```python
# Train totalVI on CITE-seq data
totalvi_model = scvi.model.TOTALVI(adata)
totalvi_model.train()

# RNA differential expression
rna_de = totalvi_model.differential_expression(
    groupby="cell_type",
    group1="T cells",
    group2="B cells",
    protein_expression=False  # Default
)

# Protein differential expression
protein_de = totalvi_model.differential_expression(
    groupby="cell_type",
    group1="T cells",
    group2="B cells",
    protein_expression=True
)

print(f"DE genes: {rna_de['is_de_fdr_0.05'].sum()}")
print(f"DE proteins: {protein_de['is_de_fdr_0.05'].sum()}")
```

### 差异化可访问性 (PeakVI)

```python
# Train PeakVI on ATAC-seq data
peakvi_model = scvi.model.PEAKVI(atac_adata)
peakvi_model.train()

# Differential accessibility
da = peakvi_model.differential_accessibility(
    groupby="cell_type",
    group1="T cells",
    group2="B cells"
)

# Same interpretation as DE
```

## 处理特殊情况

### 低细胞计数组

```python
# Increase posterior samples for stability
de = model.differential_expression(
    groupby="cell_type",
    group1="rare_type",  # e.g., 50 cells
    group2="common_type",  # e.g., 5000 cells
    n_samples=10000
)
```

### 不平衡的比较

```python
# When groups have very different sizes
# Use change mode to avoid tiny effects

de = model.differential_expression(
    groupby="condition",
    group1="rare_condition",
    group2="common_condition",
    mode="change",
    delta=0.5
)
```

### 多重测试修正

```python
# Already included via FDP control
# But can apply additional corrections

from statsmodels.stats.multitest import multipletests

# Bonferroni correction (very conservative)
_, pvals_corrected, _, _ = multipletests(
    1 / (de["bayes_factor"] + 1),
    method="bonferroni"
)
```

## 性能考虑因素

### 速度优化

```python
# Faster DE testing for large datasets
de = model.differential_expression(
    groupby="cell_type",
    group1="T cells",
    n_samples=1000,  # Reduce samples
    batch_size=512    # Increase batch size
)
```

### 内存管理

```python
# For very large datasets
# Test one comparison at a time rather than all pairwise

cell_types = adata.obs["cell_type"].unique()
for ct in cell_types:
    de = model.differential_expression(
        groupby="cell_type",
        group1=ct
    )
    # Save results
    de.to_csv(f"de_results_{ct}.csv")
```

## 最佳实践

1. **使用“改变”模式**：获得生物学上可解释的结果
2. **设置适当的增量**：基于生物学意义
3. **检查表达水平**：过滤低表达基因
4. **验证结果**：检查标记基因是否健全
5. **可视化结果**：始终绘制顶级 DE 基因
6. **报告参数**：使用的文档模式、增量、FDR
7. **考虑批量效果**：使用batch_ Correction=True
8. **多重比较**：注意测试多个组
9. **样本量**：确保每组有足够的细胞（建议> 50）
10. **生物学验证**：跟进功能实验

## 示例：完整的 DE 分析工作流程

```python
import scvi
import scanpy as sc
import matplotlib.pyplot as plt

# 1. Train model
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
model = scvi.model.SCVI(adata)
model.train()

# 2. Perform DE analysis
de_results = model.differential_expression(
    groupby="cell_type",
    group1="Disease_T_cells",
    group2="Healthy_T_cells",
    mode="change",
    delta=0.5,
    fdr_target=0.05
)

# 3. Filter and analyze
sig_genes = de_results[de_results["is_de_fdr_0.05"]]
upreg = sig_genes[sig_genes["lfc_mean"] > 0].sort_values("lfc_mean", ascending=False)
downreg = sig_genes[sig_genes["lfc_mean"] < 0].sort_values("lfc_mean")

print(f"Significant genes: {len(sig_genes)}")
print(f"Upregulated: {len(upreg)}")
print(f"Downregulated: {len(downreg)}")

# 4. Visualize top genes
top_genes = upreg.head(10).index.tolist() + downreg.head(10).index.tolist()

sc.pl.violin(
    adata[adata.obs["cell_type"].isin(["Disease_T_cells", "Healthy_T_cells"])],
    keys=top_genes,
    groupby="cell_type",
    rotation=90
)

# 5. Functional enrichment (using external tools)
# E.g., g:Profiler, DAVID, or gprofiler-official Python package
upreg_genes = upreg.head(100).index.tolist()
# Perform pathway analysis...

# 6. Save results
de_results.to_csv("de_results_disease_vs_healthy.csv")
upreg.to_csv("upregulated_genes.csv")
downreg.to_csv("downregulated_genes.csv")
```