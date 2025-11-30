<!-- 此文件由机器翻译自 coverage.md -->

# Uniwig 的覆盖率分析

uniwig 模块根据测序数据生成覆盖轨迹，提供基因组区间到覆盖图谱的高效转换。

## 覆盖轨道生成

从 BED 文件创建覆盖轨道：

```python
import gtars

# Generate coverage from BED file
coverage = gtars.uniwig.coverage_from_bed("fragments.bed")

# Generate coverage with specific resolution
coverage = gtars.uniwig.coverage_from_bed("fragments.bed", resolution=10)

# Generate strand-specific coverage
fwd_coverage = gtars.uniwig.coverage_from_bed("fragments.bed", strand="+")
rev_coverage = gtars.uniwig.coverage_from_bed("fragments.bed", strand="-")
```

## CLI 用法

从命令行生成覆盖轨迹：

<<<代码块_1>>>

## 使用覆盖率数据

### 访问覆盖范围值

查询特定位置的覆盖率：

<<<代码块_2>>>

### 覆盖操作

对覆盖轨道执行操作：

<<<代码块_3>>>

## 输出格式

Uniwig 支持多种输出格式：

### 假发格式

标准摆动格式：
<<<代码块_4>>>

### BigWig 格式

用于高效存储和访问的二进制格式：
<<<代码块_5>>>

### BedGraph 格式

可变覆盖率的灵活格式：
<<<代码块_6>>>

## 用例

### ATAC-seq 分析

生成染色质可及性概况：

```python
# Generate ATAC-seq coverage
atac_fragments = gtars.RegionSet.from_bed("atac_fragments.bed")
coverage = gtars.uniwig.coverage_from_bed("atac_fragments.bed", resolution=1)

# Identify accessible regions
peaks = coverage.call_peaks(threshold=10)
```

### ChIP-seq 峰可视化

为 ChIP-seq 数据创建覆盖轨迹：

```bash
# Generate coverage for visualization
gtars uniwig generate --input chip_seq_fragments.bed \
                      --output chip_coverage.bw \
                      --format bigwig
```

### RNA-seq 覆盖率

计算 RNA-seq 的读取覆盖率：

```python
# Generate strand-specific RNA-seq coverage
fwd = gtars.uniwig.coverage_from_bed("rnaseq.bed", strand="+")
rev = gtars.uniwig.coverage_from_bed("rnaseq.bed", strand="-")

# Export for IGV
fwd.to_bigwig("rnaseq_fwd.bw")
rev.to_bigwig("rnaseq_rev.bw")
```

### 差异覆盖分析

比较样本之间的覆盖率：

```python
# Generate coverage for two samples
control = gtars.uniwig.coverage_from_bed("control.bed")
treatment = gtars.uniwig.coverage_from_bed("treatment.bed")

# Compute fold change
fold_change = treatment.divide(control)

# Find differential regions
diff_regions = fold_change.find_regions(threshold=2.0)
```

## 性能优化

- 针对数据规模使用适当的分辨率
- 建议用于大型数据集的 BigWig 格式
- 可用于多个染色体的并行处理
- 大文件的内存高效流式传输