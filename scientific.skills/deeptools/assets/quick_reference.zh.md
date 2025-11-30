<!-- 此文件由机器翻译自 quick_reference.md -->

# deepTools 快速参考

## 最常用的命令

### BAM 到 bigWig（标准化）
```bash
bamCoverage --bam input.bam --outFileName output.bw \
    --normalizeUsing RPGC --effectiveGenomeSize 2913022398 \
    --binSize 10 --numberOfProcessors 8
```

### 比较两个 BAM 文件
<<<代码块_1>>>

### 相关热图
<<<代码块_2>>>

### TSS 周围的热图
<<<代码块_3>>>

### ChIP 富集检查
<<<代码块_4>>>

## 有效基因组大小

|有机体 |组装|尺寸|
|----------|----------|------|
|人类 | HG38 | 2913022398 |
|鼠标|毫米10 | 2652783500 |
|飞翔| DM6 | 142573017 |

## 常用归一化方法

- **RPGC**：1×基因组覆盖率（需要-- effectiveGenomeSize）
- **CPM**：每百万计数（对于固定垃圾箱）
- **RPKM**：每百万每 kb 的读取数（对于基因）

## 典型工作流程

1. **QC**：绘图指纹、绘图相关性
2. **覆盖率**：标准化后的 bamCoverage
3. **比较**：治疗组与对照组的 bamCompare
4. **可视化**：computeMatrix→plotHeatmap/plotProfile