<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 生物蟒蛇
描述：“用于分子生物学的主要 Python 工具包。首选基于 Python 的 PubMed/NCBI 查询 (Bio.Entrez)、序列操作、文件解析（FASTA、GenBank、FASTQ、PDB）、高级 BLAST 工作流程、结构、系统发育。要快速 BLAST，请使用 gget。对于直接 REST API，请使用 pubmed-database。”
---

# Biopython：Python 中的计算分子生物学

## 概述

Biopython 是一套全面的免费 Python 生物计算工具。它提供序列操作、文件 I/O、数据库访问、结构生物信息学、系统发育学和许多其他生物信息学任务的功能。当前版本是**Biopython 1.85**（2025年1月发布），支持Python 3并需要NumPy。

## 何时使用此技能

在以下情况下使用此技能：

- 处理生物序列（DNA、RNA 或蛋白质）
- 读取、写入或转换生物文件格式（FASTA、GenBank、FASTQ、PDB、mmCIF 等）
- 通过 Entrez 访问 NCBI 数据库（GenBank、PubMed、Protein、Gene 等）
- 运行 BLAST 搜索或解析 BLAST 结果
- 执行序列比对（成对或多序列比对）
- 从 PDB 文件分析蛋白质结构
- 创建、操作或可视化系统发育树
- 寻找序列基序或分析基序模式
- 计算序列统计数据（GC含量、分子量、解链温度等）
- 执行结构生物信息学任务
- 处理群体遗传学数据
- 任何其他计算分子生物学任务

## 核心能力

Biopython 被组织成模块化子包，每个子包都针对特定的生物信息学领域：

1. **序列处理** - Bio.Seq 和 Bio.SeqIO 用于序列操作和文件 I/O
2. **比对分析** - Bio.Align 和 Bio.AlignIO 用于成对和多序列比对
3. **数据库访问** - Bio.Entrez 用于以编程方式访问 NCBI 数据库
4. **BLAST 操作** - Bio.Blast 用于运行和解析 BLAST 搜索
5. **结构生物信息学** - 用于处理 3D 蛋白质结构的 Bio.PDB
6. **系统发育** - Bio.Phylo 用于系统发育树操作和可视化
7. **高级功能** - 基序、群体遗传学、序列实用程序等

## 安装和设置

使用 pip 安装 Biopython（需要 Python 3 和 NumPy）：

```python
uv pip install biopython
```

对于 NCBI 数据库访问，请始终设置您的电子邮件地址（NCBI 要求）：

<<<代码块_1>>>

## 使用此技能

该技能提供了按功能区域组织的全面文档。执行任务时，请查阅相关参考文档：

### 1. 序列处理（Bio.Seq & Bio.SeqIO）

**参考：** `references/sequence_io.md`

用于：
- 创建和操纵生物序列
- 读写序列文件（FASTA、GenBank、FASTQ等）
- 文件格式之间的转换
- 从大文件中提取序列
- 序列翻译、转录和反向互补
- 使用 SeqRecord 对象

**简单示例：**
<<<代码块_2>>>

### 2. 比对分析（Bio.Align & Bio.AlignIO）

**参考：** `references/alignment.md`

用于：
- 成对序列比对（全局和局部）
- 读取和写入多个序列比对
- 使用替换矩阵（BLOSUM、PAM）
- 计算对齐统计数据
- 自定义对齐参数

**简单示例：**
<<<代码块_3>>>

### 3. 数据库访问 (Bio.Entrez)

**参考：** `references/databases.md`

用于：
- 搜索 NCBI 数据库（PubMed、GenBank、Protein、Gene 等）
- 下载序列和记录
- 获取出版物信息
- 跨数据库查找相关记录
- 具有适当速率限制的批量下载

**简单示例：**
<<<代码块_4>>>

### 4. BLAST 操作 (Bio.Blast)

**参考：** `references/blast.md`

用于：
- 通过 NCBI Web 服务运行 BLAST 搜索
- 运行本地 BLAST 搜索
- 解析 BLAST XML 输出
- 按E值或身份过滤结果
- 提取命中序列

**简单示例：**
<<<代码块_5>>>

### 5.结构生物信息学（Bio.PDB）

**参考：** `references/structure.md`

用于：
- 解析PDB和mmCIF结构文件
- 浏览蛋白质结构层次结构（SMCRA：结构/模型/链/残基/原子）
- 计算距离、角度和二面角
- 二级结构分配（DSSP）
- 结构叠加和RMSD计算
- 从结构中提取序列

**简单示例：**
<<<代码块_6>>>

### 6.系统发育学（Bio.Phylo）

**参考：** `references/phylogenetics.md`

用于：
- 读写系统发育树（Newick、NEXUS、phyloXML）
- 从距离矩阵或路线构建树
- 树木操作（修剪、重新生根、梯形化）
- 计算系统发育距离
- 创建共识树
- 可视化树木

**简单示例：**
```python
from Bio import Phylo

# Read and visualize tree
tree = Phylo.read("tree.nwk", "newick")
Phylo.draw_ascii(tree)

# Calculate distance
distance = tree.distance("Species_A", "Species_B")
print(f"Distance: {distance:.3f}")
```

### 7. 高级功能

**参考：** `references/advanced.md`

用于：
- **序列基序** (Bio.motifs) - 查找和分析基序模式
- **群体遗传学** (Bio.PopGen) - GenePop 文件、Fst 计算、Hardy-Weinberg 测试
- **序列实用程序** (Bio.SeqUtils) - GC 含量、熔解温度、分子量、蛋白质分析
- **限制性分析** (Bio.Restriction) - 寻找限制性内切酶位点
- **聚类** (Bio.Cluster) - K-means 和层次聚类
- **基因组图** (GenomeDiagram) - 可视化基因组特征

**简单示例：**
```python
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.Seq import Seq

seq = Seq("ATCGATCGATCG")
print(f"GC content: {gc_fraction(seq):.2%}")
print(f"Molecular weight: {molecular_weight(seq, seq_type='DNA'):.2f} g/mol")
```

## 一般工作流程指南

### 阅读文档

当用户询问特定的 Biopython 任务时：

1. **根据任务描述识别相关模块**
2. **使用读取工具读取相应的参考文件**
3. **提取相关代码模式**并使其适应用户的特定需求
4. **在任务需要时组合多个模块**

参考文件的搜索模式示例：
```bash
# Find information about specific functions
grep -n "SeqIO.parse" references/sequence_io.md

# Find examples of specific tasks
grep -n "BLAST" references/blast.md

# Find information about specific concepts
grep -n "alignment" references/alignment.md
```

### 编写 Biopython 代码

编写 Biopython 代码时遵循以下原则：

1. **显式导入模块**
   ```python
   from Bio import SeqIO, Entrez
   from Bio.Seq import Seq
   ```

2. **使用 NCBI 数据库时设置 Entrez 电子邮件**
   ```python
   Entrez.email = "your.email@example.com"
   ```

3. **使用适当的文件格式** - 检查哪种格式最适合任务
   ```python
   # Common formats: "fasta", "genbank", "fastq", "clustal", "phylip"
   ```

4. **正确处理文件** - 使用后关闭句柄或使用上下文管理器
   ```python
   with open("file.fasta") as handle:
       records = SeqIO.parse(handle, "fasta")
   ```

5. **对大文件使用迭代器** - 避免将所有内容加载到内存中
   ```python
   for record in SeqIO.parse("large_file.fasta", "fasta"):
       # Process one record at a time
   ```

6. **优雅地处理错误** - 网络操作和文件解析可能会失败
   ```python
   try:
       handle = Entrez.efetch(db="nucleotide", id=accession)
   except HTTPError as e:
       print(f"Error: {e}")
   ```

## 常见模式

### 模式 1：从 GenBank 获取序列

```python
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"

# Fetch sequence
handle = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

print(f"Description: {record.description}")
print(f"Sequence length: {len(record.seq)}")
```

### 模式 2：序列分析管道

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

for record in SeqIO.parse("sequences.fasta", "fasta"):
    # Calculate statistics
    gc = gc_fraction(record.seq)
    length = len(record.seq)

    # Find ORFs, translate, etc.
    protein = record.seq.translate()

    print(f"{record.id}: {length} bp, GC={gc:.2%}")
```

### 模式 3：BLAST 并获取热门歌曲

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO

Entrez.email = "your.email@example.com"

# Run BLAST
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
blast_record = NCBIXML.read(result_handle)

# Get top hit accessions
accessions = [aln.accession for aln in blast_record.alignments[:5]]

# Fetch sequences
for acc in accessions:
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    print(f">{record.description}")
```

### 模式 4：根据序列构建系统发育树

```python
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Read alignment
alignment = AlignIO.read("alignment.fasta", "fasta")

# Calculate distances
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

# Build tree
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

# Visualize
Phylo.draw_ascii(tree)
```

## 最佳实践

1. **在编写代码之前务必阅读相关参考文档**
2. **使用grep搜索参考文件**特定函数或示例
3. **解析前验证文件格式**
4. **优雅地处理丢失的数据** - 并非所有记录都具有所有字段
5. **缓存下载的数据** - 不要重复下载相同的序列
6. **遵守 NCBI 速率限制** - 使用 API 密钥和适当的延迟
7. **在处理大文件之前使用小数据集进行测试**
8. **保持 Biopython 更新**以获得最新功能和错误修复
9. **使用适当的遗传密码表**进行翻译
10. **文件分析参数**可重复性

## 常见问题故障排除

### 问题：“找不到记录器‘Bio.Entrez’的处理程序”
**解决方案：** 这只是一个警告。设置 Entrez.email 来抑制它。

### 问题：来自 NCBI 的“HTTP 错误 400”
**解决方案：** 检查 ID/登录是否有效且格式正确。

### 问题：解析文件时出现“ValueError: EOF”
**解决方案：** 验证文件格式与指定的格式字符串匹配。

### 问题：对齐失败并显示“序列长度不同”
**解决方案：** 在使用 AlignIO 或 MultipleSeqAlignment 之前确保序列已对齐。

### 问题：BLAST 搜索速度很慢
**解决方案：** 使用本地BLAST进行大规模搜索，或者缓存结果。

### 问题：PDB 解析器警告
**解决方案：** 使用 `PDBParser(QUIET=True)` 抑制警告，或调查结构质量。

## 其他资源

- **官方文档**：https://biopython.org/docs/latest/
- **教程**：https://biopython.org/docs/latest/Tutorial/
- **食谱**：https://biopython.org/docs/latest/Tutorial/（高级示例）
- **GitHub**：https://github.com/biopython/biopython
- **邮件列表**：biopython@biopython.org

## 快速参考

要在参考文件中查找信息，请使用以下搜索模式：

```bash
# Search for specific functions
grep -n "function_name" references/*.md

# Find examples of specific tasks
grep -n "example" references/sequence_io.md

# Find all occurrences of a module
grep -n "Bio.Seq" references/*.md
```

## 总结

Biopython 为计算分子生物学提供了全面的工具。使用该技能时：

1. **识别任务域**（序列、比对、数据库、BLAST、结构、系统发育或高级）
2. **查阅 `references/` 目录中相应的参考文件**
3. **根据具体用例调整代码示例**
4. **在复杂工作流程需要时组合多个模块**
5. **遵循文件处理、错误检查和数据管理的最佳实践**

模块化参考文档确保每个主要 Biopython 功能都有详细的、可搜索的信息。