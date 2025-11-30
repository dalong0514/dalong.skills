<!-- 此文件由机器翻译自 verified-workflows.md -->

# 已验证的工作流程

## 概述
Latch 验证工作流程是由 Latch 工程师开发和维护的生产就绪型预构建生物信息学管道。这些工作流程被顶级制药公司和生物技术公司用于研究和发现。

## 在 Python SDK 中可用

`latch.verified` 模块提供从 Python 代码对经过验证的工作流程的编程访问。

### 导入已验证的工作流程

```python
from latch.verified import (
    bulk_rnaseq,
    deseq2,
    mafft,
    trim_galore,
    alphafold,
    colabfold
)
```

## 核心验证工作流程

### 批量 RNA 序列分析

**对齐和量化：**
<<<代码块_1>>>

**特点：**
- 使用 FastQC 读取质量控制
- 适配器修剪
- 与 STAR 或 HISAT2 对齐
- 使用 featureCounts 进行基因级定量
- MultiQC报告生成

### 差异表达分析

**DESeq2：**
<<<代码块_2>>>

**特点：**
- 标准化和方差稳定
- 差异表达测试
- MA地块和火山地块
- PCA可视化
- 带注释的结果表

### 路径分析

**富集分析：**
<<<代码块_3>>>

**支持的数据库：**
- 基因本体论（GO）
- KEGG途径
- 反应组
- 维基路径
- MSigDB 集合

### 序列比对

**MAFFT 多序列比对：**
<<<代码块_4>>>

**特点：**
- 多种对齐算法（FFT-NS-1、FFT-NS-2、G-INS-i、L-INS-i）
- 自动算法选择
- 支持大对齐
- 多种输出格式

### 适配器和质量修整

**修剪丰富：**
<<<代码块_5>>>

**特点：**
- 自动适配器检测
- 优质修剪
- FastQC集成
- 支持单端和双端

## 蛋白质结构预测

### 阿尔法折叠

**标准阿尔法折叠：**
<<<代码块_6>>>

**特点：**
- 单体和多聚体预测
- 基于模板的建模选项
- MSA生成
- 置信度指标（pLDDT、PAE）
- PDB结构输出

**模型预设：**
- `monomer`：单蛋白质链
- `monomer_casp14`：CASP14竞赛版本
- `monomer_ptm`：具有 pTM 信心
- `multimer`：蛋白质复合物

### ColabFold

**优化的 AlphaFold 替代方案：**
```python
from latch.verified import colabfold

structure = colabfold(
    sequence_fasta=LatchFile("latch:///data/protein.fasta"),
    num_models=5,
    use_amber_relax=True,
    output_dir="latch:///results/colabfold"
)
```

**特点：**
- 比标准 AlphaFold 更快
- 基于MMseqs2的MSA生成
- 多模型预测
- 琥珀色的放松
- 按置信度排名

**优点：**
- MSA 生成速度加快 3-5 倍
- 降低计算成本
- 与 AlphaFold 相似的精度

## 单细胞分析

### ArchR (scATAC-seq)

**染色质可及性分析：**
```python
from latch.verified import archr

results = archr(
    fragments_file=LatchFile("latch:///data/fragments.tsv.gz"),
    genome="hg38",
    output_dir="latch:///results/archr"
)
```

**特点：**
- 箭头文件生成
- 质量控制指标
- 降维
- 聚类
- 高峰通话
- 主题丰富

### scVelo（RNA 速度）

**RNA 速度分析：**
```python
from latch.verified import scvelo

results = scvelo(
    adata_file=LatchFile("latch:///data/adata.h5ad"),
    mode="dynamical",
    output_dir="latch:///results/scvelo"
)
```

**特点：**
- 剪接/未剪接定量
- 速度估计
- 动态建模
- 轨迹推断
- 可视化

###emptyDropsR（单元调用）

**空滴检测：**
```python
from latch.verified import emptydrops

filtered_matrix = emptydrops(
    raw_matrix_dir=LatchDir("latch:///data/raw_feature_bc_matrix"),
    fdr_threshold=0.01
)
```

**特点：**
- 区分细胞和空液滴
- 基于 FDR 的阈值
- 环境 RNA 去除
- 兼容10X数据

## 基因编辑分析

### CRISPResso2

**CRISPR 编辑评估：**
```python
from latch.verified import crispresso2

results = crispresso2(
    fastq_r1=LatchFile("latch:///data/sample_R1.fastq.gz"),
    amplicon_sequence="AGCTAGCTAG...",
    guide_rna="GCTAGCTAGC",
    output_dir="latch:///results/crispresso"
)
```

**特点：**
- 插入缺失定量
- 碱基编辑分析
- 主要编辑分析
- HDR量化
- 等位基因频率图

## 系统发育学

### 系统发育树构建

```python
from latch.verified import phylogenetics

tree = phylogenetics(
    alignment_file=LatchFile("latch:///data/aligned.fasta"),
    method="maximum_likelihood",
    bootstrap_replicates=1000,
    output_dir="latch:///results/phylo"
)
```

**特点：**
- 多种建树方法
- 引导程序支持
- 树可视化
- 型号选择

## 工作流程集成

### 在自定义管道中使用经过验证的工作流程

```python
from latch import workflow, small_task
from latch.verified import bulk_rnaseq, deseq2
from latch.types import LatchFile, LatchDir

@workflow
def complete_rnaseq_analysis(
    fastq_files: List[LatchFile],
    metadata: LatchFile,
    output_dir: LatchDir
) -> LatchFile:
    """
    Complete RNA-seq analysis pipeline using verified workflows
    """
    # Run alignment for each sample
    aligned_samples = []
    for fastq in fastq_files:
        result = bulk_rnaseq(
            fastq_r1=fastq,
            reference_genome="hg38",
            output_dir=output_dir
        )
        aligned_samples.append(result)

    # Aggregate counts and run differential expression
    count_matrix = aggregate_counts(aligned_samples)
    deseq_results = deseq2(
        count_matrix=count_matrix,
        sample_metadata=metadata,
        design_formula="~ condition"
    )

    return deseq_results
```

## 最佳实践

### 何时使用经过验证的工作流程

**将经过验证的工作流程用于：**
1. 标准分析管道
2. 成熟的方法
3. 生产就绪分析
4. 可重复的研究
5.经过验证的生物信息学工具

**构建自定义工作流程：**
1. 新颖的分析方法
2. 自定义预处理步骤
3. 与专有工具集成
4. 实验管线
5.高度专业化的工作流程

### 结合验证和自定义

```python
from latch import workflow, small_task
from latch.verified import alphafold
from latch.types import LatchFile

@small_task
def preprocess_sequence(raw_fasta: LatchFile) -> LatchFile:
    """Custom preprocessing"""
    # Custom logic here
    return processed_fasta

@small_task
def postprocess_structure(pdb_file: LatchFile) -> LatchFile:
    """Custom post-analysis"""
    # Custom analysis here
    return analysis_results

@workflow
def custom_structure_pipeline(input_fasta: LatchFile) -> LatchFile:
    """
    Combine custom steps with verified AlphaFold
    """
    # Custom preprocessing
    processed = preprocess_sequence(raw_fasta=input_fasta)

    # Use verified AlphaFold
    structure = alphafold(
        sequence_fasta=processed,
        model_preset="monomer_ptm"
    )

    # Custom post-processing
    results = postprocess_structure(pdb_file=structure)

    return results
```

## 访问工作流程文档

### 平台内文档

每个经过验证的工作流程包括：
- 参数说明
- 输入/输出规格
- 方法详情
- 引文信息
- 用法示例

### 查看可用的工作流程

```python
from latch.verified import list_workflows

# List all available verified workflows
workflows = list_workflows()

for workflow in workflows:
    print(f"{workflow.name}: {workflow.description}")
```

## 版本管理

### 工作流程版本

已验证的工作流程已进行版本控制和维护：
- 错误修复和改进
- 添加了新功能
- 保持向后兼容性
- 版本固定可用
### 使用特定版本

```python
from latch.verified import bulk_rnaseq

# Use specific version
results = bulk_rnaseq(
    fastq_r1=input_file,
    reference_genome="hg38",
    workflow_version="2.1.0"
)
```

## 支持和更新

### 获取帮助

- **文档**：https://docs.latch.bio
- **Slack 社区**：Latch SDK 工作区
- **支持**：support@latch.bio
- **GitHub 问题**：报告错误并请求功能

### 工作流程更新

经过验证的工作流程会定期更新：
- 工具版本升级
- 性能改进
- 错误修复
- 新功能

订阅更新通知的发行说明。

## 常见用例

### 完整的 RNA-seq 研究

```python
# 1. Quality control and alignment
aligned = bulk_rnaseq(fastq=samples)

# 2. Differential expression
deg = deseq2(counts=aligned)

# 3. Pathway enrichment
pathways = pathway_enrichment(genes=deg)
```

### 蛋白质结构分析

```python
# 1. Predict structure
structure = alphafold(sequence=protein_seq)

# 2. Custom analysis
results = analyze_structure(pdb=structure)
```

### 单细胞工作流程

```python
# 1. Filter cells
filtered = emptydrops(matrix=raw_counts)

# 2. RNA velocity
velocity = scvelo(adata=filtered)
```