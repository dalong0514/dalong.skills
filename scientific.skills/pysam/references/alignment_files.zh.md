<!-- 此文件由机器翻译自 alignment_files.md -->

# 使用对齐文件 (SAM/BAM/CRAM)

## 概述

Pysam 提供 `AlignmentFile` 类来读写包含对齐序列数据的 SAM/BAM/CRAM 格式文件。 BAM/CRAM 文件通过索引支持压缩和随机访问。

## 打开对齐文件

通过模式限定符指定格式：
- `"rb"` - 读取 BAM（二进制）
- `"r"` - 读取 SAM（文本）
- `"rc"` - 读取 CRAM（压缩）
- `"wb"` - 写入 BAM
- `"w"` - 写入 SAM
- `"wc"` - 写入 CRAM

```python
import pysam

# Reading
samfile = pysam.AlignmentFile("example.bam", "rb")

# Writing (requires template or header)
outfile = pysam.AlignmentFile("output.bam", "wb", template=samfile)
```

### 流处理

使用 `"-"` 作为 stdin/stdout 操作的文件名：

<<<代码块_1>>>

**重要提示：** Pysam 不支持从真正的 Python 文件对象读取/写入 - 仅支持 stdin/stdout 流。

## 对齐文件属性

**标题信息：**
- `references` - 染色体/重叠群名称列表
- `lengths` - 每个引用的相应长度
- `header` - 完整的标题作为字典

<<<代码块_2>>>

## 阅读阅读

### fetch() - 基于区域的检索

使用**基于 0 的坐标**检索重叠指定基因组区域的读数。

<<<代码块_3>>>

**重要说明：**
- 需要索引（.bai/.crai）进行随机访问
- 返回与区域**重叠**的读数（可能超出边界）
- 使用 `until_eof=True` 进行非索引文件或顺序读取
- 默认情况下，仅返回映射读取
- 对于未映射的读取，请使用 `fetch("*")` 或 `until_eof=True`

### 多个迭代器

在同一文件上使用多个迭代器时：

<<<代码块_4>>>

如果没有 `multiple_iterators=True`，新的 fetch() 调用将重新定位文件指针并中断现有迭代器。

### count() - 计算区域中的读取次数

<<<代码块_5>>>

### count_coverage() - 每个碱基的覆盖率

返回四个数组（A、C、G、T）以及每个碱基的覆盖范围：

<<<代码块_6>>>

## AlignedSegment 对象

每次读取都表示为具有以下关键属性的 `AlignedSegment` 对象：

### 阅读信息
- `query_name` - 读取名称/ID
- `query_sequence` - 读取序列（碱基）
- `query_qualities` - 碱基质量分数（ASCII 编码）
- `query_length` - 读取的长度

### 映射信息
- `reference_name` - 染色体/重叠群名称
- `reference_start` - 起始位置（从 0 开始，包含）
- `reference_end` - 结束位置（从 0 开始，不包括）
- `mapping_quality` - MAPQ 分数
- `cigarstring` - CIGAR 字符串（例如“100M”）
- `cigartuples` - CIGAR 作为（操作、长度）元组列表

**重要提示：** `cigartuples` 格式与 SAM 规范不同。运算都是整数：
- 0 = M（匹配/不匹配）
- 1 = I（插入）
- 2 = D（删除）
- 3 = N（跳过参考）
- 4 = S（软削波）
- 5 = H（硬削波）
- 6 = P（填充）
- 7 = =（序列匹配）
- 8 = X（序列不匹配）

### 标志和状态
- `flag` - SAM 标志为整数
- `is_paired` - 读取是否配对？
- `is_proper_pair` - 是否以正确的配对方式读取？
- `is_unmapped` - 读取是否未映射？
- `mate_is_unmapped` - 伴侣是否未映射？
- `is_reverse` - 是否在反向链上读取？
- `mate_is_reverse` - 配偶是否位于反向链上？
- `is_read1` - 这是 read1 吗？
- `is_read2` - 这是read2吗？
- `is_secondary` - 是否是二次对齐？
- `is_qcfail` - 读取 QC 失败吗？
- `is_duplicate` - 读取是否重复？
- `is_supplementary` - 是否是补充对齐？

### 标签和可选字段
- `get_tag(tag)` - 获取可选字段的值
- `set_tag(tag, value)` - 设置可选字段
- `has_tag(tag)` - 检查标签是否存在
- `get_tags()` - 获取所有标签作为元组列表

```python
for read in samfile.fetch("chr1", 1000, 2000):
    if read.has_tag("NM"):
        edit_distance = read.get_tag("NM")
        print(f"{read.query_name}: NM={edit_distance}")
```

## 写入对齐文件

### 创建标题

```python
header = {
    'HD': {'VN': '1.0'},
    'SQ': [
        {'LN': 1575, 'SN': 'chr1'},
        {'LN': 1584, 'SN': 'chr2'}
    ]
}

outfile = pysam.AlignmentFile("output.bam", "wb", header=header)
```

### 创建AlignedSegment对象

```python
# Create new read
a = pysam.AlignedSegment()
a.query_name = "read001"
a.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
a.flag = 0
a.reference_id = 0  # Index into header['SQ']
a.reference_start = 100
a.mapping_quality = 20
a.cigar = [(0, 35)]  # 35M
a.query_qualities = pysam.qualitystring_to_array("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII")

# Write to file
outfile.write(a)
```

### 格式之间的转换

```python
# BAM to SAM
infile = pysam.AlignmentFile("input.bam", "rb")
outfile = pysam.AlignmentFile("output.sam", "w", template=infile)
for read in infile:
    outfile.write(read)
infile.close()
outfile.close()
```

## 堆积分析

`pileup()` 方法提供跨区域的**列**（逐位置）分析：

```python
for pileupcolumn in samfile.pileup("chr1", 1000, 2000):
    print(f"Position {pileupcolumn.pos}: coverage = {pileupcolumn.nsegments}")

    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            # Query position is the position in the read
            base = pileupread.alignment.query_sequence[pileupread.query_position]
            print(f"  {pileupread.alignment.query_name}: {base}")
```

**关键属性：**
- `pileupcolumn.pos` - 从 0 开始的参考位置
- `pileupcolumn.nsegments` - 覆盖位置的读取数
- `pileupread.alignment` - AlignedSegment 对象
- `pileupread.query_position` - 读取中的位置（无删除）
- `pileupread.is_del` - 这是删除吗？
- `pileupread.is_refskip` - 这是引用跳过（CIGAR 中的 N）吗？
**重要：** 保持迭代器引用处于活动状态。当迭代器过早超出范围时，会出现错误“迭代器完成后访问 PileupProxy”。

## 坐标系

**关键：** Pysam 使用 **基于 0 的半开**坐标（Python 约定）：
- `reference_start` 是从 0 开始的（第一个基数是 0）
- `reference_end` 是排他的（不包含在范围内）
- 1000-2000 地区包括基地 1000-1999

**例外：** `fetch()` 和 `pileup()` 中的区域字符串遵循 samtools 约定（从 1 开始）：
```python
# These are equivalent:
samfile.fetch("chr1", 999, 2000)  # Python style: 0-based
samfile.fetch("chr1:1000-2000")   # samtools style: 1-based
```

## 索引

创建BAM索引：
```python
pysam.index("example.bam")
```

或者使用命令行界面：
```python
pysam.samtools.index("example.bam")
```

## 性能提示

1. **重复查询特定区域时使用索引访问**
2. **使用`pileup()`进行按列分析**而不是重复的获取操作
3. **使用 `fetch(until_eof=True)` 顺序读取**非索引文件
4. **除非必要，否则避免使用多个迭代器**（性能成本）
5. **使用`count()`进行简单计数**而不是手动迭代计数

## 常见陷阱

1. **部分重叠：** `fetch()` 返回与区域边界重叠的读数 - 如果需要精确边界，则实施显式过滤
2. **质量得分编辑：** 修改`query_sequence`后无法就地编辑`query_qualities`。首先创建一个副本：`quals = read.query_qualities`
3. **缺少索引：** `fetch()` 没有 `until_eof=True` 需要索引文件
4. **线程安全：**虽然pysam在I/O期间释放GIL，但全面的线程安全性尚未得到充分验证
5. **迭代器范围：** 保持堆积迭代器引用处于活动状态，以避免“迭代器完成后访问 PileupProxy”错误