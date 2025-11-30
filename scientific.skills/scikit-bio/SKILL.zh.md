<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：scikit-bio
描述：“生物数据工具包。序列分析、比对、系统发育树、多样性指标（alpha/beta、UniFrac）、排序（PCoA）、PERMANOVA、FASTA/Newick I/O，用于微生物组分析。”
---

# scikit-bio

## 概述

scikit-bio 是一个用于处理生物数据的综合 Python 库。将此技能应用于涵盖序列操作、比对、系统发育学、微生物生态学和多元统计的生物信息学分析。

## 何时使用此技能

当用户执行以下操作时应使用此技能：
- 适用于生物序列（DNA、RNA、蛋白质）
- 需要读/写生物文件格式（FASTA、FASTQ、GenBank、Newick、BIOM 等）
- 执行序列比对或搜索基序
- 构建或分析系统发育树
- 计算多样性指标（alpha/beta 多样性、UniFrac 距离）
- 执行排序分析（PCoA、CCA、RDA）
- 对生物/生态数据进行统计测试（PERMANOVA、ANOSIM、Mantel）
- 分析微生物组或群落生态数据
- 使用语言模型中的蛋白质嵌入
- 需要操作生物数据表

## 核心能力

### 1. 序列操作

使用 DNA、RNA 和蛋白质数据的专用类处理生物序列。

**关键操作：**
- 读取/写入 FASTA、FASTQ、GenBank、EMBL 格式的序列
- 序列切片、串联和搜索
- 反向互补、转录（DNA→RNA）和翻译（RNA→蛋白质）
- 使用正则表达式查找图案和模式
- 计算距离（Hamming，基于 k-mer）
- 处理序列质量分数和元数据

**常见模式：**
```python
import skbio

# Read sequences from file
seq = skbio.DNA.read('input.fasta')

# Sequence operations
rc = seq.reverse_complement()
rna = seq.transcribe()
protein = rna.translate()

# Find motifs
motif_positions = seq.find_with_regex('ATG[ACGT]{3}')

# Check for properties
has_degens = seq.has_degenerates()
seq_no_gaps = seq.degap()
```

**重要说明：**
- 使用 `DNA`、`RNA`、`Protein` 类进行语法序列验证
- 对没有字母限制的通用序列使用 `Sequence` 类
- 质量分数自动从 FASTQ 文件加载到位置元数据中
- 元数据类型：序列级（ID、描述）、位置（每个碱基）、区间（区域/特征）

### 2. 序列比对

使用动态编程算法执行成对和多序列比对。

**关键能力：**
- 全局对齐（Needleman-Wunsch 与半全局变体）
- 局部对齐（Smith-Waterman）
- 可配置的评分方案（匹配/不匹配、差距罚分、替换矩阵）
- CIGAR字符串转换
- 使用`TabularMSA`进行多序列比对存储和操作

**常见模式：**
<<<代码块_1>>>

**重要说明：**
- 使用 `local_pairwise_align_ssw` 进行局部对齐（更快，基于 SSW）
- 使用 `StripedSmithWaterman` 进行蛋白质比对
- 建议对生物序列进行仿射间隙惩罚
- 可以在 scikit-bio、BioPython 和 Biotite 比对格式之间进行转换

### 3.系统发育树

构建、操作和分析代表进化关系的系统发育树。

**关键能力：**
- 从距离矩阵构建树（UPGMA、WPGMA、邻接、GME、BME）
- 树操作（修剪、重新生根、遍历）
- 距离计算（patristic、cophenetic、Robinson-Foulds）
- ASCII 可视化
- Newick 格式 I/O

**常见模式：**
<<<代码块_2>>>

**重要说明：**
- 使用 `nj()` 进行邻居连接（经典系统发育方法）
- 使用 `upgma()` 进行 UPGMA（假设分子钟）
- GME 和 BME 对于大型树具有高度可扩展性
- 树可以是有根的，也可以是无根的；某些指标需要特定的 root

### 4.多样性分析

计算微生物生态学和群落分析的 alpha 和 beta 多样性指标。

**关键能力：**
- Alpha多样性：丰富度、香农熵、辛普森指数、Faith's PD、Pielou's均匀度
- Beta 多样性：Bray-Curtis、Jaccard、加权/未加权 UniFrac、欧几里德距离
- 系统发育多样性指标（需要树输入）
- 稀疏和二次采样
- 与排序和统计测试集成

**常见模式：**
<<<代码块_3>>>

**重要说明：**
- 计数必须是代表丰度的整数，而不是相对频率
- 系统发育指标（Faith's PD、UniFrac）需要树和 OTU ID 映射
- 使用 `partial_beta_diversity()` 仅计算特定样本对
- Alpha 多样性返回 Series，beta 多样性返回 DistanceMatrix

### 5.排序方法

将高维生物数据减少到可视化的低维空间。

**关键能力：**
- 来自距离矩阵的 PCoA（主坐标分析）
- 列联表的CA（对应分析）
- 具有环境限制的CCA（规范对应分析）
- 用于线性关系的 RDA（冗余分析）
- 用于特征解释的双图投影

**常见模式：**
<<<代码块_4>>>

**重要说明：**
- PCoA 适用于任何距离/相异矩阵
- CCA 揭示了社区组成的环境驱动因素
- 排序结果包括特征值、解释的比例和样本/特征坐标
- 结果与绘图库集成（matplotlib、seaborn、plotly）

### 6. 统计测试

执行针对生态和生物数据的假设检验。

**关键能力：**
- PERMANOVA：使用距离矩阵测试组差异
- ANOSIM：群体差异的替代测试
- PERMDISP：测试群体色散的均匀性
- Mantel检验：距离矩阵之间的相关性
- Bioenv：查找与距离相关的环境变量

**常见模式：**
<<<代码块_5>>>

**重要说明：**
- 排列检验提供非参数显着性检验
- 使用 999+ 排列以获得稳健的 p 值
- PERMANOVA 对色散差异敏感；与 PERMDISP 配对
- Mantel 检验评估矩阵相关性（例如，地理距离与遗传距离）

### 7. 文件I/O和格式转换

通过自动格式检测读取和写入 19 多种生物文件格式。

**支持的格式：**
- 序列：FASTA、FASTQ、GenBank、EMBL、QSeq
- 阵营：Clustal、PHYLIP、斯德哥尔摩
- 树：纽维克
- 表：BIOM（HDF5 和 JSON）
- 距离：定界方阵
- 分析：BLAST+6/7、GFF3、排序结果
- 元数据：带验证的 TSV/CSV

**常见模式：**
<<<代码块_6>>>

**重要说明：**
- 使用大文件生成器以避免内存问题
- 指定`into`参数时可以自动检测格式
- 某些对象可以写入多种格式
- 支持带有 `verify=False` 的 stdin/stdout 管道

### 8.距离矩阵

使用统计方法创建和操作距离/相异矩阵。

**关键能力：**
- 存储对称（DistanceMatrix）或非对称（DissimilarityMatrix）数据
- 基于ID的索引和切片
- 与多样性、排序和统计测试相结合
- 读/写分隔文本格式

**常见模式：**
```python
from skbio import DistanceMatrix
import numpy as np

# Create from array
data = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
dm = DistanceMatrix(data, ids=['A', 'B', 'C'])

# Access distances
dist_ab = dm['A', 'B']
row_a = dm['A']

# Read from file
dm = DistanceMatrix.read('distances.txt')

# Use in downstream analyses
pcoa_results = pcoa(dm)
permanova_results = permanova(dm, grouping)
```

**重要说明：**
- DistanceMatrix 强制对称性和零对角线
- DissimilarityMatrix 允许不对称值
- ID 能够与元数据和生物知识集成
- 与 pandas、numpy 和 scikit-learn 兼容

### 9. 生物表

使用微生物组研究中常见的特征表（OTU/ASV 表）。

**关键能力：**
- BIOM 格式 I/O（HDF5 和 JSON）
- 与 pandas、polars、AnnData、numpy 集成
- 数据增强技术（phylomix、mixup、组合方法）
- 样本/特征过滤和标准化
- 元数据集成

**常见模式：**
```python
from skbio import Table

# Read BIOM table
table = Table.read('table.biom')

# Access data
sample_ids = table.ids(axis='sample')
feature_ids = table.ids(axis='observation')
counts = table.matrix_data

# Filter
filtered = table.filter(sample_ids_to_keep, axis='sample')

# Convert to/from pandas
df = table.to_dataframe()
table = Table.from_dataframe(df)
```

**重要说明：**
- BIOM 表是 QIIME 2 工作流程中的标准配置
- 行通常代表样本，列代表特征（OTU/ASV）
- 支持稀疏和密集表示
- 输出格式可配置（pandas/polars/numpy）

### 10. 蛋白质嵌入

使用蛋白质语言模型嵌入进行下游分析。

**关键能力：**
- 存储蛋白质语言模型（ESM、ProtTrans 等）的嵌入
- 将嵌入转换为距离矩阵
- 生成排序对象以进行可视化
- 导出到 numpy/pandas 以用于 ML 工作流程

**常见模式：**
```python
from skbio.embedding import ProteinEmbedding, ProteinVector

# Create embedding from array
embedding = ProteinEmbedding(embedding_array, sequence_ids)

# Convert to distance matrix for analysis
dm = embedding.to_distances(metric='euclidean')

# PCoA visualization of embedding space
pcoa_results = embedding.to_ordination(metric='euclidean', method='pcoa')

# Export for machine learning
array = embedding.to_array()
df = embedding.to_dataframe()
```

**重要说明：**
- 嵌入将蛋白质语言模型与传统生物信息学联系起来
- 与 scikit-bio 的距离/排序/统计生态系统兼容
- SequenceEmbedding 和 ProteinEmbedding 提供专门的功能
- 对于序列聚类、分类和可视化很有用

## 最佳实践

### 安装
```bash
uv pip install scikit-bio
```

### 性能考虑因素
- 对大型序列文件使用生成器以最大限度地减少内存使用
- 对于大规模系统发育树，优先选择 GME 或 BME 而不是 NJ
- Beta 多样性计算可以与 `partial_beta_diversity()` 并行化
- 对于大型表，BIOM 格式 (HDF5) 比 JSON 更高效

### 与生态系统整合
- 序列通过标准格式与 Biopython 互操作
- 表格与 pandas、polars 和 AnnData 集成
- 与 scikit-learn 兼容的距离矩阵
- 使用 matplotlib/seaborn/plotly 可视化排序结果
- 与 QIIME 2 工件（BIOM、树木、距离矩阵）无缝协作

### 常见工作流程
1. **微生物组多样性分析**：读取BIOM表→计算α/β多样性→排序（PCoA）→统计测试（PERMANOVA）
2. **系统发育分析**：读取序列→对齐→构建距离矩阵→构建树→计算系统发育距离
3. **序列处理**：读取FASTQ → 质量过滤器 → 修剪/清理 → 查找主题 → 翻译 → 写入 FASTA
4. **比较基因组学**：读取序列→配对比对→计算距离→构建树→分析进化枝

## 参考文档

有关详细的 API 信息、参数规范和高级使用示例，请参阅 `references/api_reference.md`，其中包含以下内容的综合文档：
- 所有功能的完整方法签名和参数
- 复杂工作流程的扩展代码示例
- 解决常见问题
- 性能优化技巧
- 与其他库的集成模式

## 其他资源

- 官方文档：https://scikit.bio/docs/latest/
- GitHub 存储库：https://github.com/scikit-bio/scikit-bio
- 论坛支持：https://forum.qiime2.org（scikit-bio 是 QIIME 2 生态系统的一部分）