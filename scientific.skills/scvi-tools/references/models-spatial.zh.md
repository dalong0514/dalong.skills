<!-- 此文件由机器翻译自 models-spatial.md -->

# 空间转录组模型

本文档涵盖了在 scvi 工具中分析空间解析转录组数据的模型。

## DestVI（使用变分推理对空间转录组进行反卷积）

**目的**：使用单细胞参考数据对空间转录组学进行多分辨率反卷积。

**主要特点**：
- 估计每个空间位置的细胞类型比例
- 使用单细胞 RNA-seq 参考进行反卷积
- 多分辨率方法（全局和局部模式）
- 考虑空间相关性
- 提供不确定性量化

**何时使用**：
- 解卷积 Visium 或类似的空间转录组学
- 拥有带有细胞类型标签的 scRNA-seq 参考数据
- 想要将细胞类型映射到空间位置
- 对细胞类型的空间组织感兴趣
- 需要细胞类型丰度的概率估计

**数据要求**：
- **空间数据**：Visium 或类似的基于点的测量（目标数据）
- **单细胞参考**：带有细胞类型注释的 scRNA-seq
- 两个数据集应该共享基因

**基本用法**：
```python
import scvi

# Step 1: Train scVI on single-cell reference
scvi.model.SCVI.setup_anndata(sc_adata, layer="counts")
sc_model = scvi.model.SCVI(sc_adata)
sc_model.train()

# Step 2: Setup spatial data
scvi.model.DESTVI.setup_anndata(
    spatial_adata,
    layer="counts"
)

# Step 3: Train DestVI using reference
model = scvi.model.DESTVI.from_rna_model(
    spatial_adata,
    sc_model,
    cell_type_key="cell_type"  # Cell type labels in reference
)
model.train(max_epochs=2500)

# Step 4: Get cell type proportions
proportions = model.get_proportions()
spatial_adata.obsm["proportions"] = proportions

# Step 5: Get cell type-specific expression
# Expression of genes specific to each cell type at each spot
ct_expression = model.get_scale_for_ct("T cells")
```

**关键参数**：
- `amortization`：摊销策略（“两者”、“潜在”、“比例”）
- `n_latent`：潜在维度（继承自 scVI 模型）

**输出**：
- `get_proportions()`：每个点的细胞类型比例
- `get_scale_for_ct(cell_type)`：细胞类型特定的表达模式
- `get_gamma()`：特定比例的基因表达缩放

**可视化**：
<<<代码块_1>>>

## 立体镜

**目的**：使用概率建模进行空间转录组学的细胞类型反卷积。

**主要特点**：
- 基于参考的反卷积
- 细胞类型比例的概率框架
- 适用于各种空间技术
- 处理基因选择和标准化

**何时使用**：
- 与 DestVI 类似，但方法更简单
- 参考空间数据去卷积
- 基本反卷积的更快替代方案

**基本用法**：
<<<代码块_2>>>

## 七巧板

**目的**：单细胞数据到空间位置的空间映射和整合。

**主要特点**：
- 将单个细胞映射到空间坐标
- 学习单细胞和空间数据之间的最佳传输
- 空间位置的基因插补
- 细胞类型映射

**何时使用**：
- 将细胞从 scRNA-seq 映射到空间位置
- 将未测量的基因归入空间数据中
- 了解单细胞分辨率的空间组织
- 整合 scRNA-seq 和空间转录组学

**数据要求**：
- 带注释的单细胞 RNA-seq 数据
- 空间转录组学数据
- 模式之间共享基因

**基本用法**：
<<<代码块_3>>>

**可视化**：
<<<代码块_4>>>

## gimVI（用于插补的高斯恒等式 Multivi）

**目的**：空间和单细胞数据之间的跨模态插补。

**主要特点**：
- 空间和单细胞数据的联合模型
- 估算空间数据中缺失的基因
- 启用跨数据集查询
- 学习共享表示

**何时使用**：
- 估算空间数据中未测量的基因
- 空间和单细胞数据集的联合分析
- 模式之间的映射

**基本用法**：
<<<代码块_5>>>

## scVIVA（空间变分自动编码器的变化）

**目的**：分析空间数据中的细胞-环境关系。

**主要特点**：
- 模拟细胞邻域和环境
- 识别与环境相关的基因表达
- 考虑空间相关结构
- 细胞间相互作用分析

**何时使用**：
- 了解空间环境如何影响细胞
- 识别利基特异性基因程序
- 细胞间相互作用研究
- 微环境分析

**数据要求**：
- 具有坐标的空间转录组学
- 细胞类型注释（可选）

**基本用法**：
<<<代码块_6>>>

## 解决VI

**目的**：通过分辨率感知建模解决空间转录组噪音。

**主要特点**：
- 考虑空间分辨率效应
- 空间数据去噪
- 多尺度分析
- 提高下游分析质量

**何时使用**：
- 嘈杂的空间数据
- 多种空间分辨率
- 分析前需要去噪
- 提高数据质量

**基本用法**：
```python
scvi.model.RESOLVI.setup_anndata(
    spatial_adata,
    layer="counts",
    spatial_key="spatial"
)

model = scvi.model.RESOLVI(spatial_adata)
model.train()

# Get denoised expression
denoised = model.get_denoised_expression()
```

## 空间转录组学的模型选择

### 目标VI
**选择时间**：
- 需要详细的反卷积参考
- 拥有高质量的scRNA-seq参考
- 想要多分辨率分析
- 需要不确定性量化
**最适合**：Visium、基于点的技术

### 立体镜
**选择时间**：
- 需要更简单、更快的反卷积
- 基本细胞类型比例估计
- 计算资源有限

**最适合**：快速反卷积任务

### 七巧板
**选择时间**：
- 想要单细胞分辨率映射
- 需要估算许多基因
- 对细胞定位感兴趣
- 首选最佳运输方式

**最适合**：详细的空间映射

### 吉姆VI
**选择时间**：
- 需要双向插补
- 空间和单细胞联合建模
- 跨数据集查询

**最适合**：积分和插补

### scVIVA
**选择时间**：
- 对蜂窝环境感兴趣
- 细胞间相互作用分析
- 邻里效应

**最适合**：微环境研究

### 解决VI
**选择时间**：
- 数据质量是一个问题
- 需要去噪
- 多尺度分析

**最适合**：噪声数据预处理

## 完整工作流程：使用 DestVI 进行空间反卷积

```python
import scvi
import scanpy as sc
import squidpy as sq

# ===== Part 1: Prepare single-cell reference =====
# Load and process scRNA-seq reference
sc_adata = sc.read_h5ad("reference_scrna.h5ad")

# QC and filtering
sc.pp.filter_genes(sc_adata, min_cells=10)
sc.pp.highly_variable_genes(sc_adata, n_top_genes=4000)

# Train scVI on reference
scvi.model.SCVI.setup_anndata(
    sc_adata,
    layer="counts",
    batch_key="batch"
)

sc_model = scvi.model.SCVI(sc_adata)
sc_model.train(max_epochs=400)

# ===== Part 2: Load spatial data =====
spatial_adata = sc.read_visium("path/to/visium")
spatial_adata.var_names_make_unique()

# QC spatial data
sc.pp.filter_genes(spatial_adata, min_cells=10)

# ===== Part 3: Run DestVI =====
scvi.model.DESTVI.setup_anndata(
    spatial_adata,
    layer="counts"
)

destvi_model = scvi.model.DESTVI.from_rna_model(
    spatial_adata,
    sc_model,
    cell_type_key="cell_type"
)

destvi_model.train(max_epochs=2500)

# ===== Part 4: Extract results =====
# Get proportions
proportions = destvi_model.get_proportions()
spatial_adata.obsm["proportions"] = proportions

# Add proportions to .obs for easy plotting
for i, ct in enumerate(sc_model.adata.obs["cell_type"].cat.categories):
    spatial_adata.obs[f"prop_{ct}"] = proportions[:, i]

# ===== Part 5: Visualization =====
# Plot specific cell types
cell_types = ["T cells", "B cells", "Macrophages"]

for ct in cell_types:
    sc.pl.spatial(
        spatial_adata,
        color=f"prop_{ct}",
        title=f"{ct} proportions",
        spot_size=150,
        cmap="viridis"
    )

# ===== Part 6: Spatial analysis =====
# Compute spatial neighbors
sq.gr.spatial_neighbors(spatial_adata)

# Spatial autocorrelation of cell types
for ct in cell_types:
    sq.gr.spatial_autocorr(
        spatial_adata,
        attr="obs",
        mode="moran",
        genes=[f"prop_{ct}"]
    )

# ===== Part 7: Save results =====
destvi_model.save("destvi_model")
spatial_adata.write("spatial_deconvolved.h5ad")
```

## 空间分析的最佳实践

1. **参考质量**：使用高质量、注释良好的 scRNA-seq 参考
2. **基因重叠**：确保参考和空间之间有足够的共享基因
3. **空间坐标**：在`.obsm["spatial"]`中正确注册空间坐标
4. **验证**：使用已知的标记基因来验证反卷积
5. **可视化**：始终在空间上可视化结果以检查生物学的合理性
6. **细胞类型粒度**：考虑适当的细胞类型分辨率
7. **计算资源**：空间模型可能会占用大量内存
8. **质量控制**：分析前过滤低质量点