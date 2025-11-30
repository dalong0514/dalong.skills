<!-- 此文件由机器翻译自 overlap.md -->

# 重叠检测和 IGD

重叠模块使用集成基因组数据库 (IGD) 数据结构提供基因组间隔之间的有效重叠检测。

## IGD 指数

IGD（集成基因组数据库）是一种用于快速基因组区间查询和重叠检测的专用数据结构。

### 建立 IGD 指数

从基因组区域文件创建索引：

```python
import gtars

# Build IGD index from BED file
igd = gtars.igd.build_index("regions.bed")

# Save index for reuse
igd.save("regions.igd")

# Load existing index
igd = gtars.igd.load_index("regions.igd")
```

### 查询重叠

有效地找到重叠区域：

<<<代码块_1>>>

## CLI 用法

Overlaprs 命令行工具提供重叠检测：

<<<代码块_2>>>

### IGD CLI 命令

构建和查询 IGD 索引：

<<<代码块_3>>>

## Python API

### 重叠检测

计算区域集之间的重叠：

<<<代码块_4>>>

### 重叠统计

计算重叠指标：

<<<代码块_5>>>

## 性能特点

IGD 提供高效的查询：
- **索引构建**：O(n log n)，其中 n 是区域数量
- **查询时间**：O(k + log n)，其中 k 是重叠数
- **内存效率**：基因组间隔的紧凑表示

## 用例

### 监管要素分析

识别基因组特征之间的重叠：

<<<代码块_6>>>

### 变体注释

注释具有重叠特征的变体：

```python
# Check which variants overlap with coding regions
variants = gtars.RegionSet.from_bed("variants.bed")
cds = gtars.RegionSet.from_bed("coding_sequences.bed")

coding_variants = variants.filter_overlapping(cds)
```

### 染色质状态分析

比较样本之间的染色质状态：

```python
# Find regions with consistent chromatin states
sample1 = gtars.RegionSet.from_bed("sample1_peaks.bed")
sample2 = gtars.RegionSet.from_bed("sample2_peaks.bed")

consistent_regions = sample1.overlap(sample2)
```