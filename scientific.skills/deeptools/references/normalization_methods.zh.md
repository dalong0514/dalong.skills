<!-- 此文件由机器翻译自 normalization_methods.md -->

# deepTools 标准化方法

本文档解释了 deepTools 中可用的各种标准化方法以及何时使用每种方法。

## 为什么要标准化？

标准化对于以下方面至关重要：
1. **不同测序深度的样本比较**
2. **考虑文库大小差异**
3. **使覆盖率值在实验中可解释**
4. **实现条件之间的公平比较**

如果没有标准化，即使真实的生物信号相同，具有 1 亿个读数的样本似乎比具有 5000 万个读数的样本具有更高的覆盖率。

---

## 可用的标准化方法

### 1. RPKM（每百万映射读取中每千碱基的读取数）

**公式：** `(Number of reads) / (Length of region in kb × Total mapped reads in millions)`

**何时使用：**
- 比较同一样本中的不同基因组区域
- 调整测序深度和区域长度
- RNA-seq基因表达分析

**适用于：** `bamCoverage`

**示例：**
```bash
bamCoverage --bam input.bam --outFileName output.bw \
    --normalizeUsing RPKM
```

**解释：** RPKM 为 10 意味着每百万映射读取中每千碱基特征有 10 个读取。

**优点：**
- 考虑区域长度和库大小
- 在基因组学中广泛使用和理解

**缺点：**
- 如果总 RNA 含量不同，则不适合比较样品
- 比较成分差异很大的样品时可能会产生误导

---

### 2. CPM（每百万映射读取计数）

**公式：** `(Number of reads) / (Total mapped reads in millions)`

**也称为：** RPM（每百万次读取）

**何时使用：**
- 比较不同样本的相同基因组区域
- 当区域长度恒定或不相关时
- ChIP-seq、ATAC-seq、DNase-seq 分析

**适用于：** `bamCoverage`、`bamCompare`

**示例：**
<<<代码块_1>>>

**解释：** CPM 为 5 表示该 bin 中每百万映射读取有 5 个读取。

**优点：**
- 简单直观
- 适合比较不同测序深度的样本
- 适用于比较固定大小的垃圾箱

**缺点：**
- 不考虑区域长度
- 受高丰度区域（例如 RNA-seq 中的 rRNA）影响

---

### 3. BPM（每百万个映射读取的二进制数）

**公式：** `(Number of reads in bin) / (Sum of all reads in bins in millions)`

**与 CPM 的主要区别：** 仅考虑属于分析箱内的读数，而不是所有映射的读数。

**何时使用：**
- 与 CPM 类似，但当您想要排除分析区域之外的读取时
- 比较特定的基因组区域，同时忽略背景

**适用于：** `bamCoverage`、`bamCompare`

**示例：**
<<<代码块_2>>>

**解释：** BPM 仅考虑分箱区域中的读取。

**优点：**
- 将标准化重点放在分析区域上
- 受未分析区域读取的影响较小

**缺点：**
- 不太常用，可能更难与已发布的数据进行比较

---

### 4. RPGC（每个基因组内容的读取数）

**公式：** `(Number of reads × Scaling factor) / Effective genome size`

**缩放因子：** 计算以实现 1× 基因组覆盖率（每个碱基 1 次读取）

**何时使用：**
- 希望样本之间的覆盖率值具有可比性
- 需要可解释的绝对覆盖值
- 比较总读取计数差异很大的样本
- 具有尖峰标准化背景的 ChIP-seq

**适用于：** `bamCoverage`、`bamCompare`

**需要：** `--effectiveGenomeSize` 参数

**示例：**
<<<代码块_3>>>

**解释：** 信号值近似于覆盖深度（例如，2 的值 ≈ 2× 覆盖范围）。

**优点：**
- 产生 1× 归一化覆盖率
- 可根据基因组覆盖度进行解释
- 适合比较不同测序深度的样本

**缺点：**
- 需要了解有效基因组大小
- 假设覆盖范围一致（对于具有峰的 ChIP-seq 而言并非如此）

---

### 5.无（无标准化）

**公式：** 原始读取计数

**何时使用：**
- 初步分析
- 当样本具有相同的文库大小时（罕见）
- 下游工具何时执行标准化
- 调试或质量控制

**适用于：** 所有工具（通常是默认的）

**示例：**
<<<代码块_4>>>

**解释：** 每个 bin 的原始读取计数。

**优点：**
- 没有做出任何假设
- 对于查看原始数据很有用
- 最快的计算

**缺点：**
- 无法公平比较不同测序深度的样本
- 不适合出版人物

---

### 6.SES（选择性富集统计）

**方法：** 信号提取缩放 - 用于比较 ChIP 与对照的更复杂的方法

**何时使用：**
- 使用 bamCompare 进行 ChIP-seq 分析
- 想要复杂的背景校正
- 简单 readCount 缩放的替代方案
**仅适用于：** `bamCompare`

**示例：**
<<<代码块_5>>>

**注：** SES 专为 ChIP-seq 数据而设计，对于噪声数据可能比简单的读取计数缩放效果更好。

---

### 7. readCount（读取计数缩放）

**方法：** 按样本之间总读取计数的比率进行缩放

**何时使用：**
- 默认为 `bamCompare`
- 补偿比较中的测序深度差异
- 当您相信总读取计数反映了库大小时

**适用于：** `bamCompare`

**示例：**
<<<代码块_6>>>

**工作原理：** 如果sample1有100M的reads，sample2有50M的reads，则sample2在比较之前会缩放2倍。

---

## 归一化方法选择指南

### 对于 ChIP-seq 覆盖范围

**推荐：** RPGC 或 CPM

```bash
bamCoverage --bam chip.bam --outFileName chip.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398 \
    --extendReads 200 \
    --ignoreDuplicates
```

**推理：** 考虑测序深度差异； RPGC 提供可解释的覆盖率值。

---

### 对于 ChIP-seq 比较（治疗与对照）

**推荐：** log2 比率与 readCount 或 SES 缩放

```bash
bamCompare -b1 chip.bam -b2 input.bam -o ratio.bw \
    --operation log2 \
    --scaleFactorsMethod readCount \
    --extendReads 200 \
    --ignoreDuplicates
```

**推理：** Log2比率显示富集（正）和耗尽（负）； readCount 根据深度进行调整。

---

### 对于 RNA-seq 覆盖范围

**推荐：** CPM 或 RPKM

```bash
# Strand-specific forward
bamCoverage --bam rnaseq.bam --outFileName forward.bw \
    --normalizeUsing CPM \
    --filterRNAstrand forward

# For gene-level: RPKM accounts for gene length
bamCoverage --bam rnaseq.bam --outFileName output.bw \
    --normalizeUsing RPKM
```

**推理：** 用于比较固定宽度 bin 的 CPM；基因的 RPKM（考虑长度）。

---

### 对于 ATAC-seq

**推荐：** RPGC 或 CPM

```bash
bamCoverage --bam atac_shifted.bam --outFileName atac.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398
```

**推理：** 与 ChIP-seq 类似；希望样本之间的覆盖率具有可比性。

---

### 用于样本相关性分析

**推荐：** CPM 或 RPGC

```bash
multiBamSummary bins \
    --bamfiles sample1.bam sample2.bam sample3.bam \
    -o readCounts.npz

plotCorrelation -in readCounts.npz \
    --corMethod pearson \
    --whatToShow heatmap \
    -o correlation.png
```

**注意：** `multiBamSummary` 没有显式标准化，但相关性分析对缩放具有鲁棒性。对于差异很大的库大小，请考虑首先标准化 BAM 文件或使用具有 `multiBigwigSummary` 的 CPM 标准化 bigWig 文件。

---

## 高级标准化注意事项

### 尖峰归一化

对于使用掺入对照的实验（例如，用于 ChIP-seq 的 *果蝇* 染色质掺入）：

1. 根据峰值读数计算缩放因子
2. 使用 `--scaleFactor` 参数应用自定义缩放因子

```bash
# Calculate spike-in factor (example: 0.8)
SCALE_FACTOR=0.8

bamCoverage --bam chip.bam --outFileName chip_spikenorm.bw \
    --scaleFactor ${SCALE_FACTOR} \
    --extendReads 200
```

---

### 手动缩放因子

您可以应用自定义缩放因子：

```bash
# Apply 2× scaling
bamCoverage --bam input.bam --outFileName output.bw \
    --scaleFactor 2.0
```

---

### 染色体排除

从标准化计算中排除特定染色体：

```bash
bamCoverage --bam input.bam --outFileName output.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398 \
    --ignoreForNormalization chrX chrY chrM
```

**何时使用：** 混合性别样本中的性染色体、线粒体 DNA 或具有异常覆盖范围的染色体。

---

## 常见陷阱

### 1. 对基于 bin 的数据使用 RPKM
**问题：** RPKM 考虑了区域长度，但所有 bin 的大小相同
**解决方案：** 使用 CPM 或 RPGC 代替

### 2. 比较非标准化样本
**问题：** 具有 2× 测序深度的样本似乎具有 2× 信号
**解决方案：** 比较样本时始终进行归一化

### 3.有效基因组大小错误
**问题：** 使用 hg19 基因组大小获取 hg38 数据
**解决方案：** 仔细检查基因组组装并使用正确的大小

### 4. GC 校正后忽略重复项
**问题：** 可能会引入偏差
**解决方案：** 切勿在 `correctGCBias` 之后使用 `--ignoreDuplicates`

### 5. 在没有有效基因组大小的情况下使用RPGC
**问题：** 命令失败
**解决方案：** 始终在 RPGC 中指定 `--effectiveGenomeSize`

---

## 不同比较的标准化

### 样本内比较（不同地区）
**用途：** RPKM（考虑区域长度）

### 样本间比较（相同区域）
**使用：** CPM、RPGC 或 BPM（考虑库大小）

### 治疗与对照
**使用：** bam与log2比率和readCount/SES缩放比较

### 多样本相关性
**使用：** CPM或RPGC标准化bigWig文件，然后multiBigwigSummary

---

## 快速参考表

|方法|考虑深度|长度说明 |最适合 |命令|
|--------|--------------------|---------------------|----------|---------|
| RPKM | ✓ | ✓ | RNA-seq 基因 | `--normalizeUsing RPKM` |
|每千次展示费用 | ✓ | ✗ |固定尺寸垃圾箱 | `--normalizeUsing CPM` |
|业务流程管理| ✓ | ✗ |特定地区 | `--normalizeUsing BPM` |
|角色扮演游戏 | ✓ | ✗ |可解释的报道 | `--normalizeUsing RPGC --effectiveGenomeSize X` |
|无 | ✗ | ✗ |原始数据| `--normalizeUsing None` |
|社会服务局 | ✓ | ✗ | ChIP 比较 | `bamCompare --scaleFactorsMethod SES` |
|阅读次数 | ✓ | ✗ | ChIP 比较 | `bamCompare --scaleFactorsMethod readCount` |

---

## 进一步阅读

有关标准化理论和最佳实践的更多详细信息：
- deepTools 文档：https://deeptools.readthedocs.io/
- ChIP-seq 分析的 ENCODE 指南
- RNA-seq标准化论文（DESeq2、TMM方法）