<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： pydeseq2
描述：“差异基因表达分析 (Python DESeq2)。从批量 RNA-seq 计数、Wald 测试、FDR 校正、火山/MA 图中识别 DE 基因，用于 RNA-seq 分析。”
---

#PyDESeq2

## 概述

PyDESeq2 是 DESeq2 的 Python 实现，用于使用批量 RNA-seq 数据进行差异表达分析。设计和执行从数据加载到结果解释的完整工作流程，包括单因素和多因素设计、具有多重测试校正的 Wald 测试、可选的 apeGLM 收缩以及与 pandas 和 AnnData 的集成。

## 何时使用此技能

该技能应该在以下情况下使用：
- 分析大量 RNA-seq 计数数据以了解差异表达
- 比较实验条件之间的基因表达（例如，处理与对照）
- 执行考虑批次效应或协变量的多因素设计
- 将基于 R 的 DESeq2 工作流程转换为 Python
- 将差异表达分析集成到基于Python的管道中
- 用户提及“DESeq2”、“差异表达”、“RNA-seq 分析”或“PyDESeq2”

## 快速启动工作流程

对于想要执行标准差异表达分析的用户：

```python
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# 1. Load data
counts_df = pd.read_csv("counts.csv", index_col=0).T  # Transpose to samples × genes
metadata = pd.read_csv("metadata.csv", index_col=0)

# 2. Filter low-count genes
genes_to_keep = counts_df.columns[counts_df.sum(axis=0) >= 10]
counts_df = counts_df[genes_to_keep]

# 3. Initialize and fit DESeq2
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~condition",
    refit_cooks=True
)
dds.deseq2()

# 4. Perform statistical testing
ds = DeseqStats(dds, contrast=["condition", "treated", "control"])
ds.summary()

# 5. Access results
results = ds.results_df
significant = results[results.padj < 0.05]
print(f"Found {len(significant)} significant genes")
```

## 核心工作流程步骤

### 第 1 步：数据准备

**输入要求：**
- **计数矩阵：** 具有非负整数读取计数的样本 × 基因 DataFrame
- **元数据：**样本×变量数据框与实验因素

**常见的数据加载模式：**

<<<代码块_1>>>

**数据过滤：**

<<<代码块_2>>>

### 第 2 步：设计规范

设计公式指定了如何对基因表达进行建模。

**单因素设计：**
<<<代码块_3>>>

**多因素设计：**
<<<代码块_4>>>

**设计公式指南：**
- 使用威尔金森公式符号（R型）
- 将调整变量（例如批次）放在感兴趣的主要变量之前
- 确保变量作为元数据 DataFrame 中的列存在
- 使用适当的数据类型（离散变量的分类）

### 步骤 3：DESeq2 拟合

初始化 DeseqDataSet 并运行完整的管道：

<<<代码块_5>>>

**`deseq2()` 的作用：**
1. 计算尺寸因子（标准化）
2. 适合 Genewise 分散体
3.拟合离散趋势曲线
4. 计算色散先验
5. 适合MAP分散体（收缩）
6. 适应日志折叠变化
7. 计算库克距离（异常值检测）
8. 如果检测到异常值则进行重新调整（可选）

### 步骤 4：统计测试

执行 Wald 检验来识别差异表达基因：

<<<代码块_6>>>

**对比规格：**
- 格式：`[variable, test_level, reference_level]`
- 示例：`["condition", "treated", "control"]` 测试处理与对照
- 如果`None`，则使用设计中的最后一个系数

**结果数据框列：**
- `baseMean`：样本间的平均归一化计数
- `log2FoldChange`: Log2 条件之间的折叠变化
- `lfcSE`：LFC 的标准错误
- `stat`：Wald 检验统计量
- `pvalue`：原始 p 值
- `padj`：调整后的 p 值（通过 Benjamini-Hochberg 进行 FDR 校正）

### 步骤 5：可选的 LFC 收缩

应用收缩来减少倍数变化估计中的噪声：

```python
ds.lfc_shrink()  # Applies apeGLM shrinkage
```

**何时使用 LFC 收缩：**
- 用于可视化（火山图、热图）
- 按效应大小对基因进行排名
- 当对后续实验的基因进行优先排序时

**重要提示：** 收缩仅影响 log2FoldChange 值，而不影响统计测试结果（p 值保持不变）。使用缩小的值进行可视化，但报告未缩小的 p 值以显示显着性。

### 第 6 步：结果导出

保存结果和中间对象：

```python
import pickle

# Export results as CSV
ds.results_df.to_csv("deseq2_results.csv")

# Save significant genes only
significant = ds.results_df[ds.results_df.padj < 0.05]
significant.to_csv("significant_genes.csv")

# Save DeseqDataSet for later use
with open("dds_result.pkl", "wb") as f:
    pickle.dump(dds.to_picklable_anndata(), f)
```

## 常见分析模式

### 两组比较

标准病例对照比较：

```python
dds = DeseqDataSet(counts=counts_df, metadata=metadata, design="~condition")
dds.deseq2()

ds = DeseqStats(dds, contrast=["condition", "treated", "control"])
ds.summary()

results = ds.results_df
significant = results[results.padj < 0.05]
```

### 多重比较

对照对照测试多个治疗组：

```python
dds = DeseqDataSet(counts=counts_df, metadata=metadata, design="~condition")
dds.deseq2()

treatments = ["treatment_A", "treatment_B", "treatment_C"]
all_results = {}

for treatment in treatments:
    ds = DeseqStats(dds, contrast=["condition", treatment, "control"])
    ds.summary()
    all_results[treatment] = ds.results_df

    sig_count = len(ds.results_df[ds.results_df.padj < 0.05])
    print(f"{treatment}: {sig_count} significant genes")
```

### 考虑批次效应

技术变异的控制：

```python
# Include batch in design
dds = DeseqDataSet(counts=counts_df, metadata=metadata, design="~batch + condition")
dds.deseq2()

# Test condition while controlling for batch
ds = DeseqStats(dds, contrast=["condition", "treated", "control"])
ds.summary()
```

### 连续协变量

包括连续变量，如年龄或剂量：

```python
# Ensure continuous variable is numeric
metadata["age"] = pd.to_numeric(metadata["age"])

dds = DeseqDataSet(counts=counts_df, metadata=metadata, design="~age + condition")
dds.deseq2()

ds = DeseqStats(dds, contrast=["condition", "treated", "control"])
ds.summary()
```

## 使用分析脚本

该技能包括用于标准分析的完整命令行脚本：

```bash
# Basic usage
python scripts/run_deseq2_analysis.py \
  --counts counts.csv \
  --metadata metadata.csv \
  --design "~condition" \
  --contrast condition treated control \
  --output results/

# With additional options
python scripts/run_deseq2_analysis.py \
  --counts counts.csv \
  --metadata metadata.csv \
  --design "~batch + condition" \
  --contrast condition treated control \
  --output results/ \
  --min-counts 10 \
  --alpha 0.05 \
  --n-cpus 4 \
  --plots
```

**脚本特点：**
- 自动数据加载和验证
- 基因和样本过滤
- 完整的DESeq2管道执行
- 可定制参数的统计测试
- 结果导出（CSV、pickle）
- 可选的可视化（火山图和 MA 图）
当用户需要独立的分析工具或想要批量处理多个数据集时，请参考`scripts/run_deseq2_analysis.py`。

## 结果解释

### 识别重要基因

```python
# Filter by adjusted p-value
significant = ds.results_df[ds.results_df.padj < 0.05]

# Filter by both significance and effect size
sig_and_large = ds.results_df[
    (ds.results_df.padj < 0.05) &
    (abs(ds.results_df.log2FoldChange) > 1)
]

# Separate up- and down-regulated
upregulated = significant[significant.log2FoldChange > 0]
downregulated = significant[significant.log2FoldChange < 0]

print(f"Upregulated: {len(upregulated)}")
print(f"Downregulated: {len(downregulated)}")
```

### 排名和排序

```python
# Sort by adjusted p-value
top_by_padj = ds.results_df.sort_values("padj").head(20)

# Sort by absolute fold change (use shrunk values)
ds.lfc_shrink()
ds.results_df["abs_lfc"] = abs(ds.results_df.log2FoldChange)
top_by_lfc = ds.results_df.sort_values("abs_lfc", ascending=False).head(20)

# Sort by a combined metric
ds.results_df["score"] = -np.log10(ds.results_df.padj) * abs(ds.results_df.log2FoldChange)
top_combined = ds.results_df.sort_values("score", ascending=False).head(20)
```

### 质量指标

```python
# Check normalization (size factors should be close to 1)
print("Size factors:", dds.obsm["size_factors"])

# Examine dispersion estimates
import matplotlib.pyplot as plt
plt.hist(dds.varm["dispersions"], bins=50)
plt.xlabel("Dispersion")
plt.ylabel("Frequency")
plt.title("Dispersion Distribution")
plt.show()

# Check p-value distribution (should be mostly flat with peak near 0)
plt.hist(ds.results_df.pvalue.dropna(), bins=50)
plt.xlabel("P-value")
plt.ylabel("Frequency")
plt.title("P-value Distribution")
plt.show()
```

## 可视化指南

### 火山图

可视化显着性与效应大小：

```python
import matplotlib.pyplot as plt
import numpy as np

results = ds.results_df.copy()
results["-log10(padj)"] = -np.log10(results.padj)

plt.figure(figsize=(10, 6))
significant = results.padj < 0.05

plt.scatter(
    results.loc[~significant, "log2FoldChange"],
    results.loc[~significant, "-log10(padj)"],
    alpha=0.3, s=10, c='gray', label='Not significant'
)
plt.scatter(
    results.loc[significant, "log2FoldChange"],
    results.loc[significant, "-log10(padj)"],
    alpha=0.6, s=10, c='red', label='padj < 0.05'
)

plt.axhline(-np.log10(0.05), color='blue', linestyle='--', alpha=0.5)
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(Adjusted P-value)")
plt.title("Volcano Plot")
plt.legend()
plt.savefig("volcano_plot.png", dpi=300)
```

### MA 情节

显示倍数变化与平均表达：

```python
plt.figure(figsize=(10, 6))

plt.scatter(
    np.log10(results.loc[~significant, "baseMean"] + 1),
    results.loc[~significant, "log2FoldChange"],
    alpha=0.3, s=10, c='gray'
)
plt.scatter(
    np.log10(results.loc[significant, "baseMean"] + 1),
    results.loc[significant, "log2FoldChange"],
    alpha=0.6, s=10, c='red'
)

plt.axhline(0, color='blue', linestyle='--', alpha=0.5)
plt.xlabel("Log10(Base Mean + 1)")
plt.ylabel("Log2 Fold Change")
plt.title("MA Plot")
plt.savefig("ma_plot.png", dpi=300)
```

## 常见问题故障排除

### 数据格式问题

**问题：**“计数和元数据之间的索引不匹配”

**解决方案：** 确保样本名称完全匹配
```python
print("Counts samples:", counts_df.index.tolist())
print("Metadata samples:", metadata.index.tolist())

# Take intersection if needed
common = counts_df.index.intersection(metadata.index)
counts_df = counts_df.loc[common]
metadata = metadata.loc[common]
```

**问题：**“所有基因的计数为零”

**解决办法：** 检查数据是否需要转置
```python
print(f"Counts shape: {counts_df.shape}")
# If genes > samples, transpose is needed
if counts_df.shape[1] < counts_df.shape[0]:
    counts_df = counts_df.T
```

### 设计矩阵问题

**问题：**“设计矩阵不是满秩的”

**原因：** 混杂变量（例如，一批中所有处理过的样品）

**解决方案：** 删除混杂变量或添加交互项
```python
# Check confounding
print(pd.crosstab(metadata.condition, metadata.batch))

# Either simplify design or add interaction
design = "~condition"  # Remove batch
# OR
design = "~condition + batch + condition:batch"  # Model interaction
```

### 没有显着基因

**诊断：**
```python
# Check dispersion distribution
plt.hist(dds.varm["dispersions"], bins=50)
plt.show()

# Check size factors
print(dds.obsm["size_factors"])

# Look at top genes by raw p-value
print(ds.results_df.nsmallest(20, "pvalue"))
```

**可能的原因：**
- 效应量小
- 高生物变异性
- 样本量不足
- 技术问题（批次效应、异常值）

## 参考文档

有关此面向工作流程的指南之外的全面详细信息：

- **API 参考** (`references/api_reference.md`)：PyDESeq2 类、方法和数据结构的完整文档。当需要详细参数信息或了解对象属性时使用。

- **工作流程指南** (`references/workflow_guide.md`)：涵盖完整分析工作流程、数据加载模式、多因素设计、故障排除和最佳实践的深入指南。在处理复杂的实验设计或遇到问题时使用。

当用户需要时将这些引用加载到上下文中：
- 详细的API文档：`Read references/api_reference.md`
- 全面的工作流程示例：`Read references/workflow_guide.md`
- 故障排除指南：`Read references/workflow_guide.md`（请参阅故障排除部分）

## 重要提醒

1. **数据方向很重要：** 计数矩阵通常加载为基因×样本，但需要为样本×基因。如果需要，请始终使用 `.T` 进行转置。

2. **样本过滤：** 在分析之前删除元数据缺失的样本，以避免错误。

3. **基因过滤：** 过滤低计数基因（例如，总读数< 10）以提高功效并减少计算时间。

4. **设计公式顺序：** 将调整变量放在感兴趣的变量之前（例如，`"~batch + condition"` 而不是 `"~condition + batch"`）。

5. **LFC 收缩时间：** 在统计测试后应用收缩，并且仅用于可视化/排名目的。 P 值仍然基于未缩小的估计。

6. **结果解释：** 使用 `padj < 0.05` 来确定显着性，而不是原始 p 值。本杰明尼-霍赫伯格程序控制错误发现率。

7. **对比规范：** 格式为`[variable, test_level, reference_level]`，其中test_level与reference_level进行比较。

8. **保存中间对象：** 使用pickle保存DeseqDataSet对象以供以后使用或其他分析，而无需重新运行昂贵的拟合步骤。

## 安装和要求

```bash
uv pip install pydeseq2
```

**系统要求：**
-Python 3.10-3.11
- 熊猫 1.4.3+
- numpy 1.23.0+
- scipy 1.11.0+
- scikit学习1.1.1+
- 安达 0.8.0+

**可视化可选：**
-matplotlib
- 西博恩

## 其他资源

- **官方文档：** https://pydeseq2.readthedocs.io
- **GitHub 存储库：** https://github.com/owkin/PyDESeq2
- **出版物：** Muzellec 等人。 （2023）生物信息学，DOI：10.1093/生物信息学/btad547
- **原始 DESeq2 (R)：** Love 等人。 （2014）基因组生物学，DOI：10.1186/s13059-014-0550-8