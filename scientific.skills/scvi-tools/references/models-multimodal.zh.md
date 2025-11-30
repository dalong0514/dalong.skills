<!-- 此文件由机器翻译自 models-multimodal.md -->

# 多模式和多组学集成模型

本文档涵盖了 scvi-tools 中多种数据模式联合分析的模型。

##totalVI（总变分推理）

**目的**：CITE-seq 数据的联合分析（来自同一细胞的同时 RNA 和蛋白质测量）。

**主要特点**：
- 联合模拟基因表达和蛋白质丰度
- 学习共享的低维表示
- 能够根据 RNA 数据进行蛋白质推算
- 对两种模式执行差异表达
- 处理 RNA 和蛋白质层的批次效应

**何时使用**：
- 分析 CITE-seq 或 REAP-seq 数据
- 联合RNA + 表面蛋白测量
- 估算缺失的蛋白质
- 整合蛋白质和RNA信息
- 多批次CITE-seq集成

**数据要求**：
- AnnData 的基因表达位于 `.X` 或层中
- `.obsm["protein_expression"]` 中的蛋白质测量
- 两种模式测量相同的细胞

**基本用法**：
```python
import scvi

# Setup data - specify both RNA and protein layers
scvi.model.TOTALVI.setup_anndata(
    adata,
    layer="counts",  # RNA counts
    protein_expression_obsm_key="protein_expression",  # Protein counts
    batch_key="batch"
)

# Train model
model = scvi.model.TOTALVI(adata)
model.train()

# Get joint latent representation
latent = model.get_latent_representation()

# Get normalized values for both modalities
rna_normalized = model.get_normalized_expression()
protein_normalized = model.get_normalized_expression(
    transform_batch="batch1",
    protein_expression=True
)

# Differential expression (works for both RNA and protein)
rna_de = model.differential_expression(groupby="cell_type")
protein_de = model.differential_expression(
    groupby="cell_type",
    protein_expression=True
)
```

**关键参数**：
- `n_latent`：潜在空间维度（默认值：20）
- `n_layers_encoder`：编码器层数（默认值：1）
- `n_layers_decoder`：解码器层数（默认值：1）
- `protein_dispersion`：蛋白质分散处理（“蛋白质”或“蛋白质批次”）
- `empirical_protein_background_prior`：使用蛋白质的经验背景

**高级功能**：

**蛋白质插补**：
<<<代码块_1>>>

**去噪**：
<<<代码块_2>>>

**最佳实践**：
1. 对具有环境蛋白质的数据集使用经验蛋白质背景
2. 考虑异质蛋白质数据的蛋白质特异性分散
3.使用联合潜在空间进行聚类（比单独使用RNA更好）
4. 使用已知标记验证蛋白质插补
5. 训练前检查蛋白质 QC 指标

## MultiVI（多模态变分推理）

**目的**：配对和未配对的多组学数据的整合（例如，RNA + ATAC、配对和未配对的细胞）。

**主要特点**：
- 处理配对数据（相同单元格）和未配对数据（不同单元格）
- 整合多种模式：RNA、ATAC、蛋白质等。
- 缺失模态插补
- 学习跨模式的共享表示
- 灵活的整合策略

**何时使用**：
- 10x 多组数据（配对 RNA + ATAC）
- 整合单独的 RNA-seq 和 ATAC-seq 实验
- 有些细胞具有两种模式，有些只有一种模式
- 跨模态插补任务

**数据要求**：
- AnnData 具有多种模式
- 模态指标（每个细胞都有的测量值）
- 可以处理：
  - 具有两种模式的所有细胞（完全配对）
  - 配对和未配对细胞的混合
  - 完全不配对的数据集

**基本用法**：
<<<代码块_3>>>

**关键参数**：
- `n_genes`：基因特征数量
- `n_regions`：可访问区域的数量
- `n_latent`：潜在维度（默认值：20）

**集成场景**：

**场景 1：完全配对（10x Multiome）**：
<<<代码块_4>>>

**场景 2：部分配对**：
<<<代码块_5>>>

**场景 3：完全未配对**：
<<<代码块_6>>>

**高级用例**：

**跨模态预测**：
```python
# Predict peaks from gene expression
accessibility_from_rna = model.get_accessibility_estimates(
    indices=rna_only_cells
)

# Predict genes from accessibility
expression_from_atac = model.get_normalized_expression(
    indices=atac_only_cells
)
```

**特定模态分析**：
```python
# Separate analysis per modality
rna_subset = adata[adata.obs["modality"].str.contains("RNA")]
atac_subset = adata[adata.obs["modality"].str.contains("ATAC")]
```

## MrVI（多分辨率变分推理）

**目的**：多样本分析说明样本特定和共享的变异。

**主要特点**：
- 同时分析多个样品/条件
- 将变化分解为：
  - 共享变异（样本间常见）
  - 样本特定的变化
- 启用样本级别比较
- 识别样品特定的细胞状态

**何时使用**：
- 比较多个生物样本或条件
- 识别样本特定与共享细胞状态
- 疾病与健康样本比较
- 了解样本间的异质性
- 多捐助者研究

**基本用法**：
```python
scvi.model.MRVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
    sample_key="sample"  # Critical: defines biological samples
)

model = scvi.model.MRVI(adata, n_latent=10, n_latent_sample=5)
model.train()

# Get representations
shared_latent = model.get_latent_representation()  # Shared across samples
sample_specific = model.get_sample_specific_representation()

# Sample distance matrix
sample_distances = model.get_sample_distances()
```

**关键参数**：
- `n_latent`：共享潜在空间的维数
- `n_latent_sample`：样本特定空间的维数
- `sample_key`：定义生物样品的列

**分析工作流程**：
```python
# 1. Identify shared cell types across samples
sc.pp.neighbors(adata, use_rep="X_MrVI_shared")
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="shared_clusters")

# 2. Analyze sample-specific variation
sample_repr = model.get_sample_specific_representation()

# 3. Compare samples
distances = model.get_sample_distances()

# 4. Find sample-enriched genes
de_results = model.differential_expression(
    groupby="sample",
    group1="Disease",
    group2="Healthy"
)
```

**用例**：
- **多供体研究**：将供体效应与细胞类型变异分开
- **疾病研究**：识别疾病特异性与共享生物学
- **时间序列**：将时间与稳定变化分开
- **批次 + 生物学**：理清技术和生物学差异

## TotalVI 与 MultiVI 与 MrVI：何时使用哪个？

###总VI
**用于**：CITE-seq（RNA + 蛋白质，相同细胞）
- 配对测量
- 每个功能单一模态类型
- 焦点：蛋白质插补、联合分析

### 多VI
**用途**：多种模式（RNA + ATAC 等）
- 配对、未配对或混合
- 不同的功能类型
- 焦点：跨模态整合和插补

### VI先生
**用途**：多样本RNA-seq
- 单模态（RNA）
- 多个生物样本
- 焦点：样本级变异分解

## 集成最佳实践

### 对于 CITE-seq（totalVI）
1. **质控蛋白**：去除低质量抗体
2. **背景扣除**：使用经验背景先验
3. **联合聚类**：使用联合潜在空间，而不是单独使用RNA
4. **验证**：检查两种模式中的已知标记

### 对于多组/多模式 (MultiVI)
1. **特征过滤**：独立过滤基因和峰值
2. **平衡方式**：确保各方的合理代表性
3. **模态权重**：考虑一种模态是否占主导地位
4. **插补验证**：仔细验证插补值

### 对于多样本 (MrVI)
1. **样本定义**：仔细定义生物样本
2. **样本量**：每个样本需要足够的细胞
3. **协变量处理**：正确考虑批次与样本
4. **解释**：区分技术变异和生物变异

## 完整示例：使用totalVI 进行CITE-seq 分析

```python
import scvi
import scanpy as sc

# 1. Load CITE-seq data
adata = sc.read_h5ad("cite_seq.h5ad")

# 2. QC and filtering
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.highly_variable_genes(adata, n_top_genes=4000)

# Protein QC
protein_counts = adata.obsm["protein_expression"]
# Remove low-quality proteins

# 3. Setup totalVI
scvi.model.TOTALVI.setup_anndata(
    adata,
    layer="counts",
    protein_expression_obsm_key="protein_expression",
    batch_key="batch"
)

# 4. Train
model = scvi.model.TOTALVI(adata, n_latent=20)
model.train(max_epochs=400)

# 5. Extract joint representation
latent = model.get_latent_representation()
adata.obsm["X_totalVI"] = latent

# 6. Clustering on joint space
sc.pp.neighbors(adata, use_rep="X_totalVI")
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

# 7. Differential expression for both modalities
rna_de = model.differential_expression(
    groupby="leiden",
    group1="0",
    group2="1"
)

protein_de = model.differential_expression(
    groupby="leiden",
    group1="0",
    group2="1",
    protein_expression=True
)

# 8. Save model
model.save("totalvi_model")
```