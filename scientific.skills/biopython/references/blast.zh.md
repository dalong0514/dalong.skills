<!-- 此文件由机器翻译自 blast.md -->

# 使用 Bio.Blast 进行 BLAST 操作

## 概述

Bio.Blast 提供了运行 BLAST 搜索（本地和通过 NCBI Web 服务）和解析各种格式的 BLAST 结果的工具。该模块处理提交查询和解析输出的复杂性。

## 通过 NCBI Web 服务运行 BLAST

### Bio.Blast.NCBIWWW

`qblast()` 函数将序列提交到 NCBI 的在线 BLAST 服务：

```python
from Bio.Blast import NCBIWWW
from Bio import SeqIO

# Read sequence from file
record = SeqIO.read("sequence.fasta", "fasta")

# Run BLAST search
result_handle = NCBIWWW.qblast(
    program="blastn",           # BLAST program
    database="nt",              # Database to search
    sequence=str(record.seq)    # Query sequence
)

# Save results
with open("blast_results.xml", "w") as out_file:
    out_file.write(result_handle.read())
result_handle.close()
```

### 可用的 BLAST 程序

- **blastn** - 核苷酸与核苷酸
- **blastp** - 蛋白质与蛋白质
- **blastx** - 翻译的核苷酸与蛋白质
- **tblastn** - 蛋白质与翻译核苷酸
- **tblastx** - 翻译核苷酸与翻译核苷酸

### 通用数据库

**核苷酸数据库：**
- `nt` - 所有 GenBank+EMBL+DDBJ+PDB 序列
- `refseq_rna` - RefSeq RNA 序列

**蛋白质数据库：**
- `nr` - 所有非冗余 GenBank CDS 翻译
- `refseq_protein` - RefSeq 蛋白质序列
- `pdb` - 蛋白质数据库序列
- `swissprot` - 策划 UniProtKB/Swiss-Prot

### 高级 qblast 参数

<<<代码块_1>>>

### 使用序列文件或 ID

<<<代码块_2>>>

## 解析 BLAST 结果

### Bio.Blast.NCBIXML

NCBIXML 为 BLAST XML 输出（推荐格式）提供解析器：

<<<代码块_3>>>

### 访问 BLAST 记录数据

<<<代码块_4>>>

### 过滤结果

<<<代码块_5>>>

### 多个 BLAST 结果

对于包含多个 BLAST 结果的文件（例如，来自批量搜索）：

<<<代码块_6>>>

## 运行本地 BLAST

### 先决条件

本地 BLAST 要求：
1.安装BLAST+命令行工具
2.BLAST数据库下载到本地

### 使用命令行包装器

```python
from Bio.Blast.Applications import NcbiblastnCommandline

# Setup BLAST command
blastn_cline = NcbiblastnCommandline(
    query="input.fasta",
    db="local_database",
    evalue=0.001,
    outfmt=5,                    # XML format
    out="results.xml"
)

# Run BLAST
stdout, stderr = blastn_cline()

# Parse results
from Bio.Blast import NCBIXML
with open("results.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)
```

### 可用的命令行包装器

- `NcbiblastnCommandline` - BLASTN 包装器
- `NcbiblastpCommandline` - BLASTP 包装器
- `NcbiblastxCommandline` - BLASTX 包装器
- `NcbitblastnCommandline` - TBLASTN 包装器
- `NcbitblastxCommandline` - TBLASTX 包装器

### 创建 BLAST 数据库

```python
from Bio.Blast.Applications import NcbimakeblastdbCommandline

# Create nucleotide database
makedb_cline = NcbimakeblastdbCommandline(
    input_file="sequences.fasta",
    dbtype="nucl",
    out="my_database"
)
stdout, stderr = makedb_cline()
```

## 分析 BLAST 结果

### 提取最佳热门歌曲

```python
def get_best_hits(blast_record, num_hits=10, e_value_thresh=0.001):
    """Extract best hits from BLAST record."""
    hits = []
    for alignment in blast_record.alignments[:num_hits]:
        for hsp in alignment.hsps:
            if hsp.expect < e_value_thresh:
                hits.append({
                    'title': alignment.title,
                    'accession': alignment.accession,
                    'length': alignment.length,
                    'e_value': hsp.expect,
                    'score': hsp.score,
                    'identities': hsp.identities,
                    'align_length': hsp.align_length,
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end,
                    'sbjct_start': hsp.sbjct_start,
                    'sbjct_end': hsp.sbjct_end
                })
                break  # Only take best HSP per alignment
    return hits
```

### 计算同一性百分比

```python
def calculate_percent_identity(hsp):
    """Calculate percent identity for an HSP."""
    return (hsp.identities / hsp.align_length) * 100

# Use it
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < 0.001:
            identity = calculate_percent_identity(hsp)
            print(f"{alignment.title}: {identity:.2f}% identity")
```

### 提取命中序列

```python
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"

def fetch_hit_sequences(blast_record, num_sequences=5):
    """Fetch sequences for top BLAST hits."""
    sequences = []

    for alignment in blast_record.alignments[:num_sequences]:
        accession = alignment.accession

        # Fetch sequence from GenBank
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="fasta",
            retmode="text"
        )
        record = SeqIO.read(handle, "fasta")
        handle.close()

        sequences.append(record)

    return sequences
```

## 解析其他 BLAST 格式

### 制表符分隔输出 (outfmt 6/7)

```python
# Run BLAST with tabular output
blastn_cline = NcbiblastnCommandline(
    query="input.fasta",
    db="database",
    outfmt=6,
    out="results.txt"
)

# Parse tabular results
with open("results.txt") as f:
    for line in f:
        fields = line.strip().split('\t')
        query_id = fields[0]
        subject_id = fields[1]
        percent_identity = float(fields[2])
        align_length = int(fields[3])
        e_value = float(fields[10])
        bit_score = float(fields[11])

        print(f"{query_id} -> {subject_id}: {percent_identity}% identity, E={e_value}")
```

### 自定义输出格式

```python
# Specify custom columns (outfmt 6 with custom fields)
blastn_cline = NcbiblastnCommandline(
    query="input.fasta",
    db="database",
    outfmt="6 qseqid sseqid pident length evalue bitscore qseq sseq",
    out="results.txt"
)
```

## 最佳实践

1. **使用XML格式**进行解析（outfmt 5）-最可靠和完整
2. **保存 BLAST 结果** - 不必要时不要重新运行搜索
3. **设置适当的 E 值阈值** - 默认值为 10，但 0.001-0.01 通常更好
4. **处理速率限制** - NCBI 限制请求频率
5. **使用本地BLAST**进行大规模搜索或重复查询
6. **缓存结果** - 保存解析的数据以避免重新解析
7. **检查空结果** - 优雅地处理没有命中的情况
8. **考虑替代方案** - 对于大型数据集，考虑 DIAMOND 或其他快速对齐器
9. **批量搜索** - 尽可能一起提交多个序列
10. **按身份过滤** - 仅 E 值可能不够

## 常见用例

### 基本 BLAST 搜索和解析

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# Read query sequence
record = SeqIO.read("query.fasta", "fasta")

# Run BLAST
print("Running BLAST search...")
result_handle = NCBIWWW.qblast("blastn", "nt", str(record.seq))

# Parse results
blast_record = NCBIXML.read(result_handle)

# Display top 5 hits
print(f"\nTop 5 hits for {blast_record.query}:")
for i, alignment in enumerate(blast_record.alignments[:5], 1):
    hsp = alignment.hsps[0]
    identity = (hsp.identities / hsp.align_length) * 100
    print(f"{i}. {alignment.title}")
    print(f"   E-value: {hsp.expect}, Identity: {identity:.1f}%")
```

### 查找直向同源物

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"

# Query gene sequence
query_record = SeqIO.read("gene.fasta", "fasta")

# BLAST against specific organism
result_handle = NCBIWWW.qblast(
    "blastn",
    "nt",
    str(query_record.seq),
    entrez_query="Mus musculus[Organism]"  # Restrict to mouse
)

blast_record = NCBIXML.read(result_handle)

# Find best hit
if blast_record.alignments:
    best_hit = blast_record.alignments[0]
    print(f"Potential ortholog: {best_hit.title}")
    print(f"Accession: {best_hit.accession}")
```

### 批量 BLAST 多个序列

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# Read multiple sequences
sequences = list(SeqIO.parse("queries.fasta", "fasta"))

# Create batch results file
with open("batch_results.xml", "w") as out_file:
    for seq_record in sequences:
        print(f"Searching for {seq_record.id}...")

        result_handle = NCBIWWW.qblast("blastn", "nt", str(seq_record.seq))
        out_file.write(result_handle.read())
        result_handle.close()

# Parse batch results
with open("batch_results.xml") as result_handle:
    for blast_record in NCBIXML.parse(result_handle):
        print(f"\n{blast_record.query}: {len(blast_record.alignments)} hits")
```

### 互惠精选

```python
def reciprocal_best_hit(seq1_id, seq2_id, database="nr", program="blastp"):
    """Check if two sequences are reciprocal best hits."""
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio import Entrez

    Entrez.email = "your.email@example.com"

    # Forward BLAST
    result1 = NCBIWWW.qblast(program, database, seq1_id)
    record1 = NCBIXML.read(result1)
    best_hit1 = record1.alignments[0].accession if record1.alignments else None

    # Reverse BLAST
    result2 = NCBIWWW.qblast(program, database, seq2_id)
    record2 = NCBIXML.read(result2)
    best_hit2 = record2.alignments[0].accession if record2.alignments else None

    # Check reciprocity
    return best_hit1 == seq2_id and best_hit2 == seq1_id
```

## 错误处理

```python
from Bio.Blast import NCBIWWW, NCBIXML
from urllib.error import HTTPError

try:
    result_handle = NCBIWWW.qblast("blastn", "nt", "ATCGATCGATCG")
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()
except HTTPError as e:
    print(f"HTTP Error: {e.code}")
except Exception as e:
    print(f"Error running BLAST: {e}")
```