<!-- 此文件由机器翻译自 alignment.md -->

# 使用 Bio.Align 和 Bio.AlignIO 进行序列比对

## 概述

Bio.Align 提供使用各种算法进行成对序列比对的工具，而 Bio.AlignIO 则处理读取和写入各种格式的多个序列比对文件。

## 使用 Bio.Align 进行配对比对

### PairwiseAligner 类

`PairwiseAligner` 类使用 Needleman-Wunsch（全局）、Smith-Waterman（局部）、Gotoh（三态）和 Waterman-Smith-Beyer 算法执行成对序列比对。根据差距分数参数自动选择适当的算法。

### 创建对齐器

```python
from Bio import Align

# Create aligner with default parameters
aligner = Align.PairwiseAligner()

# Default scores (as of Biopython 1.85+):
# - Match score: +1.0
# - Mismatch score: 0.0
# - All gap scores: -1.0
```

### 自定义对齐参数

<<<代码块_1>>>

### 对齐模式

<<<代码块_2>>>

### 执行对齐

<<<代码块_3>>>

### 使用替换矩阵

<<<代码块_4>>>

### 可用的替换矩阵

常见的矩阵包括：
- **BLOSUM** 系列（BLOSUM45、BLOSUM50、BLOSUM62、BLOSUM80、BLOSUM90）
- **PAM** 系列（PAM30、PAM70、PAM250）
- **MATCH** - 简单匹配/不匹配矩阵

<<<代码块_5>>>

## 使用 Bio.AlignIO 进行多序列比对

### 阅读对齐

Bio.AlignIO 提供与 Bio.SeqIO 类似的 API，但用于比对文件：

<<<代码块_6>>>

### 支持的对齐格式

常见的格式包括：
- **clustal** - Clustal 格式
- **phylip** - PHYLIP 格式
- **phylip-relaxed** - 放松的 PHYLIP（较长的名称）
- **斯德哥尔摩** - 斯德哥尔摩格式
- **fasta** - FASTA 格式（对齐）
- **nexus** - NEXUS 格式
- **浮雕** - EMBOSS 对齐格式
- **msf** - MSF 格式
- **maf** - 多重对齐格式

### 写作对齐

```python
# Write alignment to file
AlignIO.write(alignment, "output.aln", "clustal")

# Convert between formats
count = AlignIO.convert("input.aln", "clustal", "output.phy", "phylip")
```

### 使用对齐对象

```python
from Bio import AlignIO

alignment = AlignIO.read("alignment.aln", "clustal")

# Get alignment properties
print(f"Number of sequences: {len(alignment)}")
print(f"Alignment length: {alignment.get_alignment_length()}")

# Access individual sequences
for record in alignment:
    print(f"{record.id}: {record.seq}")

# Get alignment column
column = alignment[:, 0]  # First column

# Get alignment slice
sub_alignment = alignment[:, 10:20]  # Positions 10-20

# Get specific sequence
seq_record = alignment[0]  # First sequence
```

### 对齐分析

```python
# Calculate alignment statistics
from Bio.Align import AlignInfo

summary = AlignInfo.SummaryInfo(alignment)

# Get consensus sequence
consensus = summary.gap_consensus(threshold=0.7)

# Position-specific scoring matrix (PSSM)
pssm = summary.pos_specific_score_matrix(consensus)

# Calculate information content
from Bio import motifs
motif = motifs.create([record.seq for record in alignment])
information = motif.counts.information_content()
```

## 以编程方式创建路线

### 来自 SeqRecord 对象

```python
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Create records
records = [
    SeqRecord(Seq("ACTGCTAGCTAG"), id="seq1"),
    SeqRecord(Seq("ACT-CTAGCTAG"), id="seq2"),
    SeqRecord(Seq("ACTGCTA-CTAG"), id="seq3"),
]

# Create alignment
alignment = MultipleSeqAlignment(records)
```

### 将序列添加到比对中

```python
# Start with empty alignment
alignment = MultipleSeqAlignment([])

# Add sequences (must have same length)
alignment.append(SeqRecord(Seq("ACTG"), id="seq1"))
alignment.append(SeqRecord(Seq("ACTG"), id="seq2"))

# Extend with another alignment
alignment.extend(other_alignment)
```

## 高级对齐操作

### 消除差距

```python
# Remove all gap-only columns
from Bio.Align import AlignInfo

no_gaps = []
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    if set(column) != {'-'}:  # Not all gaps
        no_gaps.append(column)
```

### 对齐排序

```python
# Sort by sequence ID
sorted_alignment = sorted(alignment, key=lambda x: x.id)
alignment = MultipleSeqAlignment(sorted_alignment)
```

### 计算成对恒等式

```python
def pairwise_identity(seq1, seq2):
    """Calculate percent identity between two sequences."""
    matches = sum(a == b for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    length = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
    return matches / length if length > 0 else 0

# Calculate all pairwise identities
for i, record1 in enumerate(alignment):
    for record2 in alignment[i+1:]:
        identity = pairwise_identity(record1.seq, record2.seq)
        print(f"{record1.id} vs {record2.id}: {identity:.2%}")
```

## 运行外部对齐工具

### Clustal Omega（通过命令行）

```python
from Bio.Align.Applications import ClustalOmegaCommandline

# Setup command
clustal_cmd = ClustalOmegaCommandline(
    infile="sequences.fasta",
    outfile="alignment.aln",
    verbose=True,
    auto=True
)

# Run alignment
stdout, stderr = clustal_cmd()

# Read result
alignment = AlignIO.read("alignment.aln", "clustal")
```

### 肌肉（通过命令行）

```python
from Bio.Align.Applications import MuscleCommandline

muscle_cmd = MuscleCommandline(
    input="sequences.fasta",
    out="alignment.aln"
)
stdout, stderr = muscle_cmd()
```

## 最佳实践

1. **选择适当的评分方案** - 对蛋白质使用 BLOSUM62，对 DNA 使用自定义评分
2. **考虑比对模式** - 全局用于相似长度序列，局部用于查找保守区域
3. **仔细设置间隙罚分** - 较高的罚分会产生更少、更长的间隙
4. **使用适当的格式** - FASTA用于简单比对，Stockholm用于丰富注释
5. **验证比对质量** - 检查保守区域和同一性百分比
6. **小心处理大对齐** - 使用切片和迭代来提高内存效率
7. **保留元数据** - 通过比对操作维护SeqRecord ID和注释

## 常见用例

### 找到最佳局部对齐

```python
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

aligner = PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 2
aligner.mismatch_score = -1

seq1 = Seq("AGCTTAGCTAGCTAGC")
seq2 = Seq("CTAGCTAGC")

alignments = aligner.align(seq1, seq2)
print(alignments[0])
```

### 蛋白质序列比对

```python
from Bio.Align import PairwiseAligner, substitution_matrices

aligner = PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

protein1 = Seq("KEVLA")
protein2 = Seq("KEVLAEQP")
alignments = aligner.align(protein1, protein2)
```

### 提取保守区域

```python
from Bio import AlignIO

alignment = AlignIO.read("alignment.aln", "clustal")

# Find columns with >80% identity
conserved_positions = []
for i in range(alignment.get_alignment_length()):
    column = alignment[:, i]
    most_common = max(set(column), key=column.count)
    if column.count(most_common) / len(column) > 0.8:
        conserved_positions.append(i)

print(f"Conserved positions: {conserved_positions}")
```