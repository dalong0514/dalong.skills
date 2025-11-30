<!-- 此文件由机器翻译自 refget.md -->

# 参考序列管理

refget 模块处理参考序列检索和摘要计算，遵循 refget 协议进行序列识别。

## 引用存储

RefgetStore 管理参考序列及其摘要：

```python
import gtars

# Create RefgetStore
store = gtars.RefgetStore()

# Add sequence
store.add_sequence("chr1", sequence_data)

# Retrieve sequence
seq = store.get_sequence("chr1")

# Get sequence digest
digest = store.get_digest("chr1")
```

## 序列摘要

计算并验证序列摘要：

<<<代码块_1>>>

## 与参考基因组整合

使用标准参考基因组：

<<<代码块_2>>>

## CLI 用法

从命令行管理参考序列：

<<<代码块_3>>>

## Refget 协议合规性

refget 模块遵循 GA4GH refget 协议：

### 摘要计算

摘要使用截断为 48 字节的 SHA-512 计算：

<<<代码块_4>>>

### 序列检索

通过摘要检索序列：

<<<代码块_5>>>

## 用例

### 参考验证

验证参考基因组的完整性：

<<<代码块_6>>>

### 序列提取

提取特定基因组区域：

```python
# Extract regions of interest
store = gtars.RefgetStore.from_fasta("hg38.fa")

regions = [
    ("chr1", 1000, 2000),
    ("chr2", 5000, 6000),
    ("chr3", 10000, 11000)
]

sequences = [store.get_subsequence(c, s, e) for c, s, e in regions]
```

### 交叉参考比较

比较不同参考文献中的序列：

```python
# Load two reference versions
hg19 = gtars.RefgetStore.from_fasta("hg19.fa")
hg38 = gtars.RefgetStore.from_fasta("hg38.fa")

# Compare digests
for chrom in hg19.chromosomes:
    digest_19 = hg19.get_digest(chrom)
    digest_38 = hg38.get_digest(chrom)
    if digest_19 != digest_38:
        print(f"{chrom} differs between hg19 and hg38")
```

## 性能说明

- 按需加载序列
- 计算后缓存摘要
- 高效的子序列提取
- 对大型基因组的内存映射文件支持