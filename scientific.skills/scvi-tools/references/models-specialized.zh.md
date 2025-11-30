<!-- 此文件由机器翻译自 models-specialized.md -->

# 专门的模态模型

本文档涵盖了 scvi 工具中专门的单细胞数据模式的模型。

## MmethylVI / MmethylANVI（甲基化分析）

**目的**：分析 DNA 甲基化的单细胞亚硫酸氢盐测序 (scBS-seq) 数据。

**主要特点**：
- 以单细胞分辨率模拟甲基化模式
- 处理甲基化数据的稀疏性
- 甲基化实验的批量校正
- 用于细胞类型注释的标签转移（MmethylANVI）

**何时使用**：
- 分析 scBS-seq 或类似的甲基化数据
- 研究跨细胞类型的 DNA 甲基化模式
- 跨批次整合甲基化数据
- 基于甲基化谱的细胞类型注释

**数据要求**：
- 甲基化计数矩阵（甲基化与每个 CpG 位点的总读数）
- 格式：细胞 × CpG 位点以及甲基化比率或计数

### 甲基六（无监督）

**基本用法**：
```python
import scvi

# Setup methylation data
scvi.model.METHYLVI.setup_anndata(
    adata,
    layer="methylation_counts",  # Methylation data
    batch_key="batch"
)

model = scvi.model.METHYLVI(adata)
model.train()

# Get latent representation
latent = model.get_latent_representation()

# Get normalized methylation values
normalized_meth = model.get_normalized_methylation()
```

### MmethylANVI（半监督细胞类型）

**基本用法**：
<<<代码块_1>>>

**关键参数**：
- `n_latent`：潜在维度
- `region_factors`：模型区域特定效果

**用例**：
- 表观遗传异质性分析
- 通过甲基化鉴定细胞类型
- 与基因表达数据整合（单独分析）
- 差异甲基化分析

## CytoVI（流式细胞术和质谱流式细胞仪）

**目的**：流式细胞术和质谱流式细胞术 (CyTOF) 数据的批量校正和整合。

**主要特点**：
- 处理基于抗体的蛋白质测量
- 修正细胞计数数据中的批次效应
- 实现跨实验的集成
- 专为高维蛋白质面板设计

**何时使用**：
- 分析流式细胞术或 CyTOF 数据
- 跨批次整合细胞计数实验
- 蛋白质面板的批量校正
- 交叉研究细胞计数整合

**数据要求**：
- 蛋白质表达矩阵（细胞×蛋白质）
- 流式细胞术或 CyTOF 测量
- 批次/实验注释

**基本用法**：
<<<代码块_2>>>

**关键参数**：
- `n_latent`：潜在空间维度
- `n_layers`：网络深度

**典型工作流程**：
<<<代码块_3>>>

## SysVI（系统级集成）

**目的**：批量效应校正，重点是保留生物变异。

**主要特点**：
- 专门的批量集成方法
- 保留生物信号，同时消除技术影响
- 专为大规模集成研究而设计

**何时使用**：
- 大规模多批次集成
- 需要保留微妙的生物变异
- 跨多项研究的系统级分析

**基本用法**：
<<<代码块_4>>>

## 破译（轨迹推断）

**目的**：单细胞数据的轨迹推断和伪时间分析。

**主要特点**：
- 学习细胞轨迹和分化路径
- 伪时间估计
- 考虑轨迹结构的不确定性
- 与 scVI 嵌入兼容

**何时使用**：
- 研究细胞分化
- 时间进程或发育数据集
- 了解细胞状态转换
- 识别发展的分支点

**基本用法**：
<<<代码块_5>>>

**可视化**：
<<<代码块_6>>>

## peRegLM（峰值调节线性模型）

**目的**：将染色质可及性与基因表达联系起来以进行调控分析。

**主要特点**：
- 将 ATAC-seq 峰与基因表达联系起来
- 确定监管关系
- 使用配对的多组数据

**何时使用**：
- 多组数据（来自相同细胞的RNA + ATAC）
- 了解基因调控
- 将峰与目标基因联系起来
- 监管网络建设

**基本用法**：
```python
# Requires paired RNA + ATAC data
scvi.model.PEREGLM.setup_anndata(
    multiome_adata,
    rna_layer="counts",
    atac_layer="atac_counts"
)

model = scvi.model.PEREGLM(multiome_adata)
model.train()

# Get peak-gene links
peak_gene_links = model.get_regulatory_links()
```

## 特定于模型的最佳实践

### 甲基VI/甲基ANVI
1. **稀疏性**：甲基化数据本质上是稀疏的；模型解释了这一点
2. **CpG选择**：过滤覆盖度极低的CpG
3. **生物学解释**：考虑基因组背景（启动子、增强子）
4. **整合**：对于多组学，单独分析然后整合结果

### CytoVI
1. **蛋白质QC**：去除低质量或无信息的蛋白质
2. **补偿**：分析前确保适当的光谱补偿
3. **批量设计**：包括生物和技术复制
4. **对照**：使用对照样品来验证批次校正

### 系统VI
1. **样本大小**：专为大规模集成而设计
2. **批次定义**：仔细定义批次结构
3. **生物验证**：验证保留的生物信号
### 破译
1. **起点**：定义轨迹起始单元（如果已知）
2. **分支**：指定期望的分支数量
3. **验证**：使用已知标记来验证伪时间
4. **集成**：与 scVI 嵌入配合良好

## 与其他模型的集成

许多专业模型可以很好地组合使用：

**甲基化+表达**：
```python
# Analyze separately, then integrate
methylvi_model = scvi.model.METHYLVI(meth_adata)
scvi_model = scvi.model.SCVI(rna_adata)

# Integrate results at analysis level
# E.g., correlate methylation and expression patterns
```

**细胞计数 + CITE-seq**：
```python
# CytoVI for flow/CyTOF
cyto_model = scvi.model.CYTOVI(cyto_adata)

# totalVI for CITE-seq
cite_model = scvi.model.TOTALVI(cite_adata)

# Compare protein measurements across platforms
```

**ATAC + RNA（多组）**：
```python
# MultiVI for joint analysis
multivi_model = scvi.model.MULTIVI(multiome_adata)

# peRegLM for regulatory links
pereglm_model = scvi.model.PEREGLM(multiome_adata)
```

## 选择专业型号

### 决策树

1. **什么数据方式？**
   - 甲基化 → 甲基VI/甲基ANVI
   - 流式/CyTOF → CytoVI
   - 轨迹→破译
   - 多批次集成 → SysVI
   - 监管链接 → peRegLM

2. **你们有标签吗？**
   - 是 → MmethylANVI（甲基化）
   - 否 → MmethylVI（甲基化）

3. **你的主要目标是什么？**
   - 批量校正 → CytoVI、SysVI
   - 轨迹/伪时间→破译
   - 峰基因链接 → peRegLM
   - 甲基化模式 → 甲基VI/ANVI

## 示例：完整的甲基化分析

```python
import scvi
import scanpy as sc

# 1. Load methylation data
meth_adata = sc.read_h5ad("methylation_data.h5ad")

# 2. QC: filter low-coverage CpG sites
sc.pp.filter_genes(meth_adata, min_cells=10)

# 3. Setup MethylVI
scvi.model.METHYLVI.setup_anndata(
    meth_adata,
    layer="methylation",
    batch_key="batch"
)

# 4. Train model
model = scvi.model.METHYLVI(meth_adata, n_latent=15)
model.train(max_epochs=400)

# 5. Get latent representation
latent = model.get_latent_representation()
meth_adata.obsm["X_MethylVI"] = latent

# 6. Clustering
sc.pp.neighbors(meth_adata, use_rep="X_MethylVI")
sc.tl.umap(meth_adata)
sc.tl.leiden(meth_adata)

# 7. Differential methylation
dm_results = model.differential_methylation(
    groupby="leiden",
    group1="0",
    group2="1"
)

# 8. Save
model.save("methylvi_model")
meth_adata.write("methylation_analyzed.h5ad")
```

## 外部工具集成

一些专用模型可作为外部包提供：

**SOLO**（双峰检测）：
```python
from scvi.external import SOLO

solo = SOLO.from_scvi_model(scvi_model)
solo.train()
doublets = solo.predict()
```

**scArches**（参考映射）：
```python
from scvi.external import SCARCHES

# For transfer learning and query-to-reference mapping
```

这些外部工具针对特定用例扩展了 scvi-tools 功能。

## 汇总表

|型号|数据类型|主要用途 |监督？ |
|--------|---------|-------------|-------------|
|甲基VI |甲基化 |无监督分析 |没有 |
|甲基ANVI |甲基化 |细胞类型注释|半|
|细胞VI |细胞计数|批量修正|没有 |
|系统VI |单链RNA测序|大规模整合|没有 |
|破译|单链RNA测序|轨迹推断 |没有 |
| peRegLM |多组学|峰基因链接|没有 |
|独奏|单链RNA测序|双峰检测|半|