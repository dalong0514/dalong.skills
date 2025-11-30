<!-- 此文件由机器翻译自 workflow_guide.md -->

# PyDESeq2 工作流程指南

本文档提供常见 PyDESeq2 分析模式的详细分步工作流程。

## 目录
1. [完全差异表达分析](#complete- Differential-express-analysis)
2. [数据加载和准备](#data-loading-and-preparation)
3. [单因素分析](#single-factor-analysis)
4. [多因素分析](#multi-factor-analysis)
5. [结果导出和可视化](#result-export-and-visualization)
6. [常见模式和最佳实践](#common-patterns-and-best-practices)
7. [疑难解答](#疑难解答)

---

## 完整的差异表达分析

### 概述
标准 PyDESeq2 分析由跨两个阶段的 12 个主要步骤组成：

**阶段 1：读取计数建模（步骤 1-7）**
- 归一化和离散度估计
原木折叠变化配件
- 异常值检测

**第 2 阶段：统计分析（步骤 8-12）**
- 瓦尔德测试
- 多重测试修正
- 可选的 LFC 收缩

### 完整工作流程代码

```python
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Load data
counts_df = pd.read_csv("counts.csv", index_col=0).T  # Transpose if needed
metadata = pd.read_csv("metadata.csv", index_col=0)

# Filter low-count genes
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
counts_df = counts_df[genes_to_keep]

# Remove samples with missing metadata
samples_to_keep = ~metadata.condition.isna()
counts_df = counts_df.loc[samples_to_keep]
metadata = metadata.loc[samples_to_keep]

# Initialize DeseqDataSet
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~condition",
    refit_cooks=True
)

# Run normalization and fitting
dds.deseq2()

# Perform statistical testing
ds = DeseqStats(
    dds,
    contrast=["condition", "treated", "control"],
    alpha=0.05,
    cooks_filter=True,
    independent_filter=True
)
ds.summary()

# Optional: Apply LFC shrinkage for visualization
ds.lfc_shrink()

# Access results
results = ds.results_df
print(results.head())
```

---

## 数据加载和准备

### 加载 CSV 文件

计数数据通常采用基因×样本格式，但需要转置：

<<<代码块_1>>>

### 从其他格式加载

**来自 TSV：**
<<<代码块_2>>>

**来自保存的泡菜：**
<<<代码块_3>>>

**来自安数据：**
<<<代码块_4>>>

### 数据过滤

**过滤低计数基因：**
<<<代码块_5>>>

**过滤缺少元数据的样本：**
<<<代码块_6>>>

**按多个条件过滤：**
```python
# Keep only samples that meet all criteria
mask = (
    ~metadata.condition.isna() &
    (metadata.batch.isin(["batch1", "batch2"])) &
    (metadata.age >= 18)
)
counts_df = counts_df.loc[mask]
metadata = metadata.loc[mask]
```

### 数据验证

**检查数据结构：**
```python
print(f"Counts shape: {counts_df.shape}")  # Should be (samples, genes)
print(f"Metadata shape: {metadata.shape}")  # Should be (samples, variables)
print(f"Indices match: {all(counts_df.index == metadata.index)}")

# Check for negative values
assert (counts_df >= 0).all().all(), "Counts must be non-negative"

# Check for non-integer values
assert counts_df.applymap(lambda x: x == int(x)).all().all(), "Counts must be integers"
```

---

## 单因素分析

### 简单的两组比较

比较处理样品与对照样品：

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Design: model expression as a function of condition
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~condition"
)

dds.deseq2()

# Test treated vs control
ds = DeseqStats(
    dds,
    contrast=["condition", "treated", "control"]
)
ds.summary()

# Results
results = ds.results_df
significant = results[results.padj < 0.05]
print(f"Found {len(significant)} significant genes")
```

### 多个成对比较

比较多个组时：

```python
# Test each treatment vs control
treatments = ["treated_A", "treated_B", "treated_C"]
all_results = {}

for treatment in treatments:
    ds = DeseqStats(
        dds,
        contrast=["condition", treatment, "control"]
    )
    ds.summary()
    all_results[treatment] = ds.results_df

# Compare results across treatments
for name, results in all_results.items():
    sig = results[results.padj < 0.05]
    print(f"{name}: {len(sig)} significant genes")
```

---

## 多因素分析

### 双因素设计

测试条件时考虑批次效应：

```python
# Design includes both batch and condition
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~batch + condition"
)

dds.deseq2()

# Test condition effect while controlling for batch
ds = DeseqStats(
    dds,
    contrast=["condition", "treated", "control"]
)
ds.summary()
```

### 交互效果

测试组间治疗效果是否不同：

```python
# Design includes interaction term
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~group + condition + group:condition"
)

dds.deseq2()

# Test the interaction term
ds = DeseqStats(dds, contrast=["group:condition", ...])
ds.summary()
```

### 连续协变量

包括连续变量，例如年龄：

```python
# Ensure age is numeric in metadata
metadata["age"] = pd.to_numeric(metadata["age"])

dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~age + condition"
)

dds.deseq2()
```

---

## 结果导出和可视化

### 保存结果

**导出为 CSV：**
```python
# Save statistical results
ds.results_df.to_csv("deseq2_results.csv")

# Save significant genes only
significant = ds.results_df[ds.results_df.padj < 0.05]
significant.to_csv("significant_genes.csv")

# Save with sorted results
sorted_results = ds.results_df.sort_values("padj")
sorted_results.to_csv("sorted_results.csv")
```

**保存 DeseqDataSet：**
```python
import pickle

# Save as AnnData for later use
with open("dds_result.pkl", "wb") as f:
    pickle.dump(dds.to_picklable_anndata(), f)
```

**加载保存的结果：**
```python
# Load results
results = pd.read_csv("deseq2_results.csv", index_col=0)

# Load AnnData
with open("dds_result.pkl", "rb") as f:
    adata = pickle.load(f)
```

### 基本可视化

**火山图：**
```python
import matplotlib.pyplot as plt
import numpy as np

results = ds.results_df.copy()
results["-log10(padj)"] = -np.log10(results.padj)

# Plot
plt.figure(figsize=(10, 6))
plt.scatter(
    results.log2FoldChange,
    results["-log10(padj)"],
    alpha=0.5,
    s=10
)
plt.axhline(-np.log10(0.05), color='red', linestyle='--', label='padj=0.05')
plt.axvline(1, color='gray', linestyle='--')
plt.axvline(-1, color='gray', linestyle='--')
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(Adjusted P-value)")
plt.title("Volcano Plot")
plt.legend()
plt.savefig("volcano_plot.png", dpi=300)
```

**MA图：**
```python
plt.figure(figsize=(10, 6))
plt.scatter(
    np.log10(results.baseMean + 1),
    results.log2FoldChange,
    alpha=0.5,
    s=10,
    c=(results.padj < 0.05),
    cmap='bwr'
)
plt.xlabel("Log10(Base Mean + 1)")
plt.ylabel("Log2 Fold Change")
plt.title("MA Plot")
plt.savefig("ma_plot.png", dpi=300)
```

---

## 常见模式和最佳实践

### 1. 数据预处理清单

运行 PyDESeq2 之前：
- ✓ 确保计数为非负整数
- ✓ 验证样本×基因方向
- ✓ 检查样本名称在计数和元数据之间是否匹配
- ✓ 删除或处理缺失的元数据值
- ✓ 过滤低计数基因（通常< 10 总读数）
- ✓ 验证实验因素是否正确编码

### 2.设计公式最佳实践

**顺序很重要：** 将调整变量放在感兴趣的变量之前
```python
# Correct: control for batch, test condition
design = "~batch + condition"

# Less ideal: condition listed first
design = "~condition + batch"
```

**对离散变量使用分类：**
```python
# Ensure proper data types
metadata["condition"] = metadata["condition"].astype("category")
metadata["batch"] = metadata["batch"].astype("category")
```

### 3. 统计测试指南

**设置适当的阿尔法：**
```python
# Standard significance threshold
ds = DeseqStats(dds, alpha=0.05)

# More stringent for exploratory analysis
ds = DeseqStats(dds, alpha=0.01)
```

**使用独立过滤：**
```python
# Recommended: filter low-power tests
ds = DeseqStats(dds, independent_filter=True)

# Only disable if you have specific reasons
ds = DeseqStats(dds, independent_filter=False)
```

### 4. LFC 收缩

**何时使用：**
- 用于可视化（火山图、热图）
- 按效应大小对基因进行排名
- 当优先考虑基因进行后续研究时

**何时不使用：**
- 用于报告统计显着性（使用未缩小的 p 值）
- 用于基因集富集分析（通常使用未缩小的值）

```python
# Save both versions
ds.results_df.to_csv("results_unshrunken.csv")
ds.lfc_shrink()
ds.results_df.to_csv("results_shrunken.csv")
```

### 5.内存管理

对于大型数据集：
```python
# Use parallel processing
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~condition",
    n_cpus=4  # Adjust based on available cores
)

# Process in batches if needed
# (split genes into chunks, analyze separately, combine results)
```

---

## 故障排除

### 错误：计数和元数据之间的索引不匹配

**问题：** 示例名称不匹配
```
KeyError: Sample names in counts and metadata don't match
```

**解决方案：**
```python
# Check indices
print("Counts samples:", counts_df.index.tolist())
print("Metadata samples:", metadata.index.tolist())

# Align if needed
common_samples = counts_df.index.intersection(metadata.index)
counts_df = counts_df.loc[common_samples]
metadata = metadata.loc[common_samples]
```

### 错误：所有基因的计数为零

**问题：** 数据可能需要转置
```
ValueError: All genes have zero total counts
```

**解决方案：**
```python
# Check data orientation
print(f"Counts shape: {counts_df.shape}")

# If genes > samples, likely needs transpose
if counts_df.shape[1] < counts_df.shape[0]:
    counts_df = counts_df.T
```

### 警告：许多基因被过滤掉

**问题：** 删除了太多低计数基因

**检查：**
```python
# See distribution of gene counts
print(counts_df.sum(axis=0).describe())

# Visualize
import matplotlib.pyplot as plt
plt.hist(counts_df.sum(axis=0), bins=50, log=True)
plt.xlabel("Total counts per gene")
plt.ylabel("Frequency")
plt.show()
```

**如果需要调整过滤：**
```python
# Try lower threshold
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 5]
```

### 错误：设计矩阵不是满秩的

**问题：** 混杂设计（例如，一批中所有处理过的样品）

**解决方案：**
```python
# Check design confounding
print(pd.crosstab(metadata.condition, metadata.batch))

# Either remove confounded variable or add interaction term
design = "~condition"  # Drop batch
# OR
design = "~condition + batch + condition:batch"  # Add interaction
```

### 问题：未发现重要基因

**可能的原因：**
1. 效应量小
2. 生物变异性高
3、样本量不足
4. 技术问题（批次效应、异常值）

**诊断：**
```python
# Check dispersion estimates
import matplotlib.pyplot as plt
dispersions = dds.varm["dispersions"]
plt.hist(dispersions, bins=50)
plt.xlabel("Dispersion")
plt.ylabel("Frequency")
plt.show()

# Check size factors (should be close to 1)
print("Size factors:", dds.obsm["size_factors"])

# Look at top genes even if not significant
top_genes = ds.results_df.nsmallest(20, "pvalue")
print(top_genes)
```

### 大型数据集上的内存错误

**解决方案：**
```python
# 1. Use fewer CPUs (paradoxically can help)
dds = DeseqDataSet(..., n_cpus=1)

# 2. Filter more aggressively
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 20]

# 3. Process in batches
# Split analysis by gene subsets and combine results
```