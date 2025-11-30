<!-- 此文件由机器翻译自 variant_files.md -->

# 使用变体文件 (VCF/BCF)

## 概述

Pysam 提供了 `VariantFile` 类来读写 VCF（变体调用格式）和 BCF（二进制 VCF）文件。这些文件包含有关遗传变异的信息，包括 SNP、插入缺失和结构变异。

## 打开变体文件

```python
import pysam

# Reading VCF
vcf = pysam.VariantFile("example.vcf")

# Reading BCF (binary, compressed)
bcf = pysam.VariantFile("example.bcf")

# Reading compressed VCF
vcf_gz = pysam.VariantFile("example.vcf.gz")

# Writing
outvcf = pysam.VariantFile("output.vcf", "w", header=vcf.header)
```

## 变体文件属性

**标题信息：**
- `header` - 包含元数据的完整 VCF 标头
- `header.contigs` - 重叠群/染色体字典
- `header.samples` - 样本名称列表
- `header.filters` - FILTER 定义字典
- `header.info` - INFO 字段定义字典
- `header.formats` - FORMAT 字段定义字典

<<<代码块_1>>>

## 读取变体记录

### 迭代所有变体

<<<代码块_2>>>

### 获取特定区域

VCF.gz 需要 tabix 索引 (.tbi) 或 BCF 需要索引：

<<<代码块_3>>>

**注意：** 在 `fetch()` 调用中使用 **基于 1 的坐标**以匹配 VCF 规范。

## VariantRecord 对象

每个变体都表示为一个 `VariantRecord` 对象：

### 职位信息
- `chrom` - 染色体/重叠群名称
- `pos` - 位置（从 1 开始）
- `start` - 起始位置（从 0 开始）
- `stop` - 停止位置（从 0 开始，不包括）
- `id` - 变体 ID（例如 rsID）

### 等位基因信息
- `ref` - 参考等位基因
- `alts` - 替代等位基因的元组
- `alleles` - 所有等位基因的元组（ref + alts）

### 质量和过滤
- `qual` - 质量得分（QUAL 字段）
- `filter` - 过滤器状态

### 信息字段

作为字典访问 INFO 字段：

<<<代码块_4>>>

### 基因型数据样本

通过 `samples` 字典访问示例数据：

<<<代码块_5>>>

**基因型代表：**
- `(0, 0)` - 纯合参考
- `(0, 1)` - 杂合子
- `(1, 1)` - 纯合子替代
- `(None, None)` - 缺失基因型
- 分相：`(0|1)` 与非分相：`(0/1)`

## 编写变体文件

### 创建标题

<<<代码块_6>>>

### 创建变体记录

```python
# Create new variant
record = outvcf.new_record()
record.chrom = "chr1"
record.pos = 100000
record.id = "rs123456"
record.ref = "A"
record.alts = ("G",)
record.qual = 30
record.filter.add("PASS")

# Set INFO fields
record.info["DP"] = 100
record.info["AF"] = (0.25,)

# Set genotype data
record.samples["sample1"]["GT"] = (0, 1)
record.samples["sample1"]["DP"] = 50
record.samples["sample2"]["GT"] = (0, 0)
record.samples["sample2"]["DP"] = 50

# Write to file
outvcf.write(record)
```

## 过滤变体

### 基本过滤

```python
# Filter by quality
for variant in vcf:
    if variant.qual >= 30:
        print(f"High quality variant: {variant.chrom}:{variant.pos}")

# Filter by depth
for variant in vcf:
    if "DP" in variant.info and variant.info["DP"] >= 20:
        print(f"High depth variant: {variant.chrom}:{variant.pos}")

# Filter by allele frequency
for variant in vcf:
    if "AF" in variant.info:
        for af in variant.info["AF"]:
            if af >= 0.01:
                print(f"Common variant: {variant.chrom}:{variant.pos}")
```

### 按基因型过滤

```python
# Find variants where sample has alternate allele
for variant in vcf:
    sample = variant.samples["sample1"]
    gt = sample["GT"]

    # Check if has alternate allele
    if gt and any(allele and allele > 0 for allele in gt):
        print(f"Sample has alt allele: {variant.chrom}:{variant.pos}")

    # Check if homozygous alternate
    if gt == (1, 1):
        print(f"Homozygous alt: {variant.chrom}:{variant.pos}")
```

### 过滤字段

```python
# Check FILTER status
for variant in vcf:
    if "PASS" in variant.filter or len(variant.filter) == 0:
        print(f"Passed filters: {variant.chrom}:{variant.pos}")
    else:
        print(f"Failed: {variant.filter.keys()}")
```

## 索引 VCF 文件

为压缩的VCF创建tabix索引：

```python
# Compress and index
pysam.tabix_index("example.vcf", preset="vcf", force=True)
# Creates example.vcf.gz and example.vcf.gz.tbi
```

或者使用 bcftools 进行 BCF：

```python
pysam.bcftools.index("example.bcf")
```

## 常见工作流程

### 提取特定样本的变体

```python
invcf = pysam.VariantFile("input.vcf")
samples_to_keep = ["sample1", "sample3"]

# Create new header with subset of samples
new_header = invcf.header.copy()
new_header.samples.clear()
for sample in samples_to_keep:
    new_header.samples.add(sample)

outvcf = pysam.VariantFile("output.vcf", "w", header=new_header)

for variant in invcf:
    # Create new record
    new_record = outvcf.new_record(
        contig=variant.chrom,
        start=variant.start,
        stop=variant.stop,
        alleles=variant.alleles,
        id=variant.id,
        qual=variant.qual,
        filter=variant.filter,
        info=variant.info
    )

    # Copy genotype data for selected samples
    for sample in samples_to_keep:
        new_record.samples[sample].update(variant.samples[sample])

    outvcf.write(new_record)
```

### 计算等位基因频率

```python
vcf = pysam.VariantFile("example.vcf")

for variant in vcf:
    total_alleles = 0
    alt_alleles = 0

    for sample_name in variant.samples:
        gt = variant.samples[sample_name]["GT"]
        if gt and None not in gt:
            total_alleles += 2
            alt_alleles += sum(1 for allele in gt if allele > 0)

    if total_alleles > 0:
        af = alt_alleles / total_alleles
        print(f"{variant.chrom}:{variant.pos} AF={af:.4f}")
```

### 将 VCF 转换为汇总表

```python
import csv

vcf = pysam.VariantFile("example.vcf")

with open("variants.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "DP"])

    for variant in vcf:
        writer.writerow([
            variant.chrom,
            variant.pos,
            variant.id or ".",
            variant.ref,
            ",".join(variant.alts) if variant.alts else ".",
            variant.qual or ".",
            variant.info.get("DP", ".")
        ])
```

## 性能提示

1. **使用BCF格式**比VCF有更好的压缩效果和更快的访问速度
2. **使用tabix索引文件**以实现高效的区域查询
3. **尽早过滤**以减少对不相关变体的处理
4. **有效使用INFO字段** - 访问前检查是否存在
5. **创建VCF文件时批量写入操作**

## 常见陷阱

1. **坐标系：** VCF使用从1开始的坐标，但VariantRecord.start是从0开始的
2. **丢失数据：** 在访问之前始终检查 INFO/FORMAT 字段是否存在
3. **基因型元组：** 基因型是元组，而不是列表——处理缺失数据的 None 值
4. **等位基因索引：** 在基因型 (0, 1) 中，0=REF，1=第一个 ALT，2=第二个 ALT，等等。
5. **索引要求：** 基于区域的 `fetch()` 需要 VCF.gz 的 tabix 索引
6. **标头修改：** 对样本进行子集化时，正确更新标头并复制 FORMAT 字段