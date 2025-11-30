<!-- 此文件由机器翻译自 workflows.md -->

# deepTools 常用工作流程

本文档提供了常见 deepTools 分析的完整工作流程示例。

## ChIP-seq 质量控制工作流程

ChIP-seq 实验的完整质量控制评估。

### 第 1 步：初始相关性评估

比较重复和样品以验证实验质量：

```bash
# Generate coverage matrix across genome
multiBamSummary bins \
    --bamfiles Input1.bam Input2.bam ChIP1.bam ChIP2.bam \
    --labels Input_rep1 Input_rep2 ChIP_rep1 ChIP_rep2 \
    -o readCounts.npz \
    --numberOfProcessors 8

# Create correlation heatmap
plotCorrelation \
    -in readCounts.npz \
    --corMethod pearson \
    --whatToShow heatmap \
    --plotFile correlation_heatmap.png \
    --plotNumbers

# Generate PCA plot
plotPCA \
    -in readCounts.npz \
    -o PCA_plot.png \
    -T "PCA of ChIP-seq samples"
```

**预期结果：**
- 重复应该聚集在一起
- 输入样本应与 ChIP 样本不同

---

### 第 2 步：覆盖范围和深度评估

<<<代码块_1>>>

**解释：** 评估测序深度是否足以进行下游分析。

---

### 步骤 3：片段大小验证（双端）

<<<代码块_2>>>

**预期结果：** 片段大小应符合文库制备方案（ChIP-seq 通常为 200-600bp）。

---

### 步骤 4：GC 偏差检测和校正

<<<代码块_3>>>

**注意：** 仅在观察到显着偏差时才正确。请勿将 `--ignoreDuplicates` 与 GC 更正的文件一起使用。

---

### 步骤 5：ChIP 信号强度评估

<<<代码块_4>>>

**释义：**
- 强大的 ChIP：累积曲线急剧上升
- 弱丰富：曲线接近对角线（类似输入）

---

## ChIP-seq 分析工作流程

从 BAM 文件到出版物质量可视化的完整工作流程。

### 第 1 步：生成标准化覆盖轨道

<<<代码块_5>>>

---

### 步骤 2：创建 Log2 比率轨迹

<<<代码块_6>>>

**结果：** Log2 比率轨迹显示富集（正值）和消耗（负值）。

---

### 步骤 3：计算 TSS 周围的矩阵

```bash
# Prepare data for heatmap/profile around transcription start sites
computeMatrix reference-point \
    --referencePoint TSS \
    --scoreFileName ChIP_coverage.bw \
    --regionsFileName genes.bed \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --binSize 10 \
    --sortRegions descend \
    --sortUsing mean \
    --outFileName matrix_TSS.gz \
    --outFileNameMatrix matrix_TSS.tab \
    --numberOfProcessors 8
```

---

### 步骤 4：生成热图

```bash
# Create heatmap around TSS
plotHeatmap \
    --matrixFile matrix_TSS.gz \
    --outFileName heatmap_TSS.png \
    --colorMap RdBu \
    --whatToShow 'plot, heatmap and colorbar' \
    --zMin -3 --zMax 3 \
    --yAxisLabel "Genes" \
    --xAxisLabel "Distance from TSS (bp)" \
    --refPointLabel "TSS" \
    --heatmapHeight 15 \
    --kmeans 3
```

---

### 步骤 5：生成剖面图

```bash
# Create meta-profile around TSS
plotProfile \
    --matrixFile matrix_TSS.gz \
    --outFileName profile_TSS.png \
    --plotType lines \
    --perGroup \
    --colors blue \
    --plotTitle "ChIP-seq signal around TSS" \
    --yAxisLabel "Average signal" \
    --xAxisLabel "Distance from TSS (bp)" \
    --refPointLabel "TSS"
```

---

### 步骤 6：峰值富集

```bash
# Calculate enrichment in peak regions
plotEnrichment \
    --bamfiles Input.bam ChIP.bam \
    --BED peaks.bed \
    --labels Input ChIP \
    --plotFile enrichment.png \
    --outRawCounts enrichment_counts.tab \
    --extendReads 200 \
    --ignoreDuplicates
```

---

## RNA-seq 覆盖工作流程

为 RNA-seq 数据生成链特异性覆盖轨迹。

### 前进链

```bash
bamCoverage \
    --bam rnaseq.bam \
    --outFileName forward_coverage.bw \
    --filterRNAstrand forward \
    --normalizeUsing CPM \
    --binSize 1 \
    --numberOfProcessors 8
```

### 反向链

```bash
bamCoverage \
    --bam rnaseq.bam \
    --outFileName reverse_coverage.bw \
    --filterRNAstrand reverse \
    --normalizeUsing CPM \
    --binSize 1 \
    --numberOfProcessors 8
```

**重要提示：** 不要使用 `--extendReads` 进行 RNA-seq（会延伸到剪接点）。

---

## 多样本比较工作流程

比较多个 ChIP-seq 样本（例如，不同的条件或时间点）。

### 第 1 步：生成覆盖率文件

```bash
# For each sample
for sample in Control_ChIP Treated_ChIP; do
    bamCoverage \
        --bam ${sample}.bam \
        --outFileName ${sample}.bw \
        --normalizeUsing RPGC \
        --effectiveGenomeSize 2913022398 \
        --binSize 10 \
        --extendReads 200 \
        --ignoreDuplicates \
        --numberOfProcessors 8
done
```

---

### 步骤 2：计算多样本矩阵

```bash
computeMatrix scale-regions \
    --scoreFileName Control_ChIP.bw Treated_ChIP.bw \
    --regionsFileName genes.bed \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    --regionBodyLength 3000 \
    --binSize 10 \
    --sortRegions descend \
    --sortUsing mean \
    --outFileName matrix_multi.gz \
    --numberOfProcessors 8
```

---

### 步骤 3：多样本热图

```bash
plotHeatmap \
    --matrixFile matrix_multi.gz \
    --outFileName heatmap_comparison.png \
    --colorMap Blues \
    --whatToShow 'plot, heatmap and colorbar' \
    --samplesLabel Control Treated \
    --yAxisLabel "Genes" \
    --heatmapHeight 15 \
    --kmeans 4
```

---

### 步骤 4：多样本分析

```bash
plotProfile \
    --matrixFile matrix_multi.gz \
    --outFileName profile_comparison.png \
    --plotType lines \
    --perGroup \
    --colors blue red \
    --samplesLabel Control Treated \
    --plotTitle "ChIP-seq signal comparison" \
    --startLabel "TSS" \
    --endLabel "TES"
```

---

## ATAC-seq 工作流程

带有 Tn5 偏移校正的 ATAC-seq 数据的专业工作流程。

### 第 1 步：移位读数以进行 Tn5 校正

```bash
alignmentSieve \
    --bam atacseq.bam \
    --outFile atacseq_shifted.bam \
    --ATACshift \
    --minFragmentLength 38 \
    --maxFragmentLength 2000 \
    --ignoreDuplicates
```

---

### 第 2 步：生成覆盖轨迹

```bash
bamCoverage \
    --bam atacseq_shifted.bam \
    --outFileName atacseq_coverage.bw \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398 \
    --binSize 1 \
    --numberOfProcessors 8
```

---

### 步骤 3：片段大小分析

```bash
bamPEFragmentSize \
    --bamfiles atacseq.bam \
    --histogram fragmentSizes_atac.png \
    --maxFragmentLength 1000
```

**预期模式：** 核小体阶梯，峰值位于 ~50bp（无核小体）、~200bp（单核小体）、~400bp（双核小体）。

---

## 峰区域分析工作流程

特别分析峰区域的 ChIP-seq 信号。

### 步骤 1：峰值矩阵

```bash
computeMatrix reference-point \
    --referencePoint center \
    --scoreFileName ChIP_coverage.bw \
    --regionsFileName peaks.bed \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 2000 \
    --binSize 10 \
    --outFileName matrix_peaks.gz \
    --numberOfProcessors 8
```

---

### 步骤 2：峰值热图

```bash
plotHeatmap \
    --matrixFile matrix_peaks.gz \
    --outFileName heatmap_peaks.png \
    --colorMap YlOrRd \
    --refPointLabel "Peak Center" \
    --heatmapHeight 15 \
    --sortUsing max
```

---

## 常见问题故障排除

### 问题：内存不足
**解决方案：** 使用`--region`参数单独处理染色体：
```bash
bamCoverage --bam input.bam -o chr1.bw --region chr1
```

### 问题：BAM 索引缺失
**解决方案：** 在运行 deepTools 之前索引 BAM 文件：
```bash
samtools index input.bam
```

### 问题：处理速度慢
**解决方案：**增加`--numberOfProcessors`：
```bash
# Use 8 cores instead of default
--numberOfProcessors 8
```

### 问题：bigWig 文件太大
**解决方案：** 增加垃圾箱大小：
```bash
--binSize 50  # or larger (default is 10-50)
```

---

## 性能提示

1. **使用多个处理器：** 始终将 `--numberOfProcessors` 设置为可用内核
2. **进程区域：** 使用 `--region` 进行测试或内存有限的环境
3. **调整 bin 大小：** 较大的 bin = 更快的处理速度和更小的文件
4. **预过滤BAM文件：** 使用`alignmentSieve`创建过滤后的BAM文件一次，然后重复使用
5. **在 bedGraph 上使用 bigWig：** bigWig 格式经过压缩，处理速度更快

---

## 最佳实践

1. **始终首先检查 QC：** 在继续之前运行相关性、覆盖率和指纹分析
2. **文档参数：** 保存命令行以实现重现性
3. **使用一致的归一化：** 在比较中的样本之间应用相同的归一化方法
4. **验证参考基因组匹配：** 确保 BAM 文件和区域文件使用相同的基因组构建
5. **检查链方向：** 对于 RNA-seq，验证正确的链方向
6. **先在小区域上测试：** 使用 `--region chr1:1-1000000` 测试参数
7. **保留中间文件：** 保存矩阵以使用不同设置重新生成绘图