<!-- 此文件由机器翻译自 multiparametric.md -->

# 多参数成像

## 概述

PathML 为多参数成像技术提供专门支持，以单细胞分辨率同时测量多个标记物。这些技术包括 CODEX、Vectra 多重免疫荧光、MERFISH 和其他空间蛋白质组学和转录组学平台。 PathML 处理特定于每种技术的独特数据结构、处理要求和量化工作流程。

## 支持的技术

### CODEX（通过 indEXing 进行协同检测）
- 循环免疫荧光成像
- 同时检测 40 多个蛋白质标记
- 单细胞空间蛋白质组学
- 使用抗体条形码进行多周期采集

### 威达北极星
- 多光谱多重免疫荧光
- 每张幻灯片 6-8 个标记
- 光谱分解
- 全玻片扫描

### MERFISH（多重误差稳健 FISH）
- 空间转录组学
- 数百至数千个基因
- 单分子分辨率
- 纠错条形码

### 其他平台
- CycIF（循环免疫荧光）
- IMC（成像质量细胞仪）
- MIBI（多重离子束成像）

## CODEX 工作流程

### 加载 CODEX 数据

CODEX 数据通常组织在来自多个采集周期的多通道图像堆栈中：

```python
from pathml.core import CODEXSlide

# Load CODEX dataset
codex_slide = CODEXSlide(
    path='path/to/codex_directory',
    stain='IF',  # Immunofluorescence
    backend='bioformats'
)

# Inspect channels and cycles
print(f"Number of channels: {codex_slide.num_channels}")
print(f"Channel names: {codex_slide.channel_names}")
print(f"Number of cycles: {codex_slide.num_cycles}")
print(f"Image shape: {codex_slide.shape}")
```

**CODEX目录结构：**
<<<代码块_1>>>

### CODEX 预处理管道

CODEX 数据处理的完整管道：

<<<代码块_2>>>

### CollapseRunsCODEX

将多周期 CODEX 采集整合到单个多通道图像中：

<<<代码块_3>>>

**参数：**
- `z_slice`：从 z 堆栈中提取哪个焦平面（通常是中间切片）
- `run_order`：循环顺序；无自动检测功能
- `method`：如何组合多个周期的通道（'max'、'mean'、'median'）
- `background_subtract`：是否减去背景荧光

**输出：** 带有所有标记（H、W、C）的单个多通道图像

### 使用 Mesmer 进行细胞分割

DeepCell Mesmer 为多参数成像提供准确的细胞分割：

<<<代码块_4>>>

**选择细胞质通道：**
- **CD45**：泛白细胞标记物（适用于免疫丰富的组织）
- **panCK**：全细胞角蛋白（有益于上皮组织）
- **CD298/b2m**：通用膜标记物
- **组合**：平均多个膜标记

**隔间选项：**
- `'whole-cell'`：完整细胞分割（细胞核 + 细胞质）
- `'nuclear'`：仅核分割
- `'cytoplasm'`：仅限细胞质区室

### 远程分段

使用 DeepCell 云 API 进行分割，无需本地 GPU：

<<<代码块_5>>>

### 标记定量

从分割图像中提取单细胞标记表达：

<<<代码块_6>>>

**输出：** AnnData 对象具有：
- `adata.X`：标记表达矩阵（细胞×标记）
- `adata.obs`：单元格元数据（单元格 ID、坐标、面积等）
- `adata.var`：标记元数据
- `adata.obsm['spatial']`：单元格质心坐标

### 与 AnnData 集成

将多个 CODEX 幻灯片处理为统一的 AnnData 对象：

```python
from pathml.core import SlideDataset
import anndata as ad

# Process multiple slides
slide_paths = ['slide1', 'slide2', 'slide3']
dataset = SlideDataset(
    [CODEXSlide(p, stain='IF') for p in slide_paths]
)

# Run pipeline on all slides
dataset.run(codex_pipeline, distributed=True, n_workers=8)

# Combine into single AnnData
adatas = []
for slide in dataset:
    adata = slide.cell_data
    adata.obs['slide_id'] = slide.name
    adatas.append(adata)

# Concatenate
combined_adata = ad.concat(adatas, join='outer', label='batch', keys=slide_paths)

# Save for downstream analysis
combined_adata.write('codex_dataset.h5ad')
```

## Vectra 工作流程

### 加载 Vectra 数据

Vectra 以专有的 `.qptiff` 格式存储数据：

```python
from pathml.core import SlideData, SlideType

# Load Vectra slide
vectra_slide = SlideData.from_slide(
    'path/to/slide.qptiff',
    backend=SlideType.VectraQPTIFF
)

# Access spectral channels
print(f"Channels: {vectra_slide.channel_names}")
```

### Vectra 预处理

```python
from pathml.preprocessing import Pipeline, CollapseRunsVectra, SegmentMIF, QuantifyMIF

vectra_pipeline = Pipeline([
    # 1. Process Vectra multi-channel data
    CollapseRunsVectra(
        wavelengths=[520, 540, 570, 620, 670, 780],  # Emission wavelengths
        unmix=True,  # Apply spectral unmixing
        autofluorescence_correction=True
    ),

    # 2. Cell segmentation
    SegmentMIF(
        nuclear_channel='DAPI',
        cytoplasm_channel='FITC',
        model='mesmer',
        image_resolution=0.5
    ),

    # 3. Quantification
    QuantifyMIF(
        segmentation_mask_name='cell_segmentation',
        markers=['DAPI', 'CD3', 'CD8', 'PD1', 'PDL1', 'panCK'],
        output_format='anndata'
    )
])

vectra_pipeline.run(vectra_slide)
```

## 下游分析

### 单元格类型注释

根据标记表达注释单元格：

```python
import anndata as ad
import numpy as np

# Load quantified data
adata = ad.read_h5ad('codex_dataset.h5ad')

# Define cell types by marker thresholds
def annotate_cell_types(adata, thresholds):
    cell_types = np.full(adata.n_obs, 'Unknown', dtype=object)

    # T cells: CD3+
    cd3_pos = adata[:, 'CD3'].X.flatten() > thresholds['CD3']
    cell_types[cd3_pos] = 'T cell'

    # CD4 T cells: CD3+ CD4+ CD8-
    cd4_tcells = (
        (adata[:, 'CD3'].X.flatten() > thresholds['CD3']) &
        (adata[:, 'CD4'].X.flatten() > thresholds['CD4']) &
        (adata[:, 'CD8'].X.flatten() < thresholds['CD8'])
    )
    cell_types[cd4_tcells] = 'CD4 T cell'

    # CD8 T cells: CD3+ CD8+ CD4-
    cd8_tcells = (
        (adata[:, 'CD3'].X.flatten() > thresholds['CD3']) &
        (adata[:, 'CD8'].X.flatten() > thresholds['CD8']) &
        (adata[:, 'CD4'].X.flatten() < thresholds['CD4'])
    )
    cell_types[cd8_tcells] = 'CD8 T cell'

    # B cells: CD20+
    b_cells = adata[:, 'CD20'].X.flatten() > thresholds['CD20']
    cell_types[b_cells] = 'B cell'

    # Macrophages: CD68+
    macrophages = adata[:, 'CD68'].X.flatten() > thresholds['CD68']
    cell_types[macrophages] = 'Macrophage'

    # Tumor cells: panCK+
    tumor = adata[:, 'panCK'].X.flatten() > thresholds['panCK']
    cell_types[tumor] = 'Tumor'

    return cell_types

# Apply annotation
thresholds = {
    'CD3': 0.5,
    'CD4': 0.4,
    'CD8': 0.4,
    'CD20': 0.3,
    'CD68': 0.3,
    'panCK': 0.5
}

adata.obs['cell_type'] = annotate_cell_types(adata, thresholds)

# Visualize cell type composition
import matplotlib.pyplot as plt
cell_type_counts = adata.obs['cell_type'].value_counts()
plt.figure(figsize=(10, 6))
cell_type_counts.plot(kind='bar')
plt.xlabel('Cell Type')
plt.ylabel('Count')
plt.title('Cell Type Composition')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()
```

### 聚类

用于识别细胞群的无监督聚类：

```python
import scanpy as sc

# Preprocessing for clustering
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

# PCA
sc.tl.pca(adata, n_comps=50)

# Neighborhood graph
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

# UMAP embedding
sc.tl.umap(adata)

# Leiden clustering
sc.tl.leiden(adata, resolution=0.5)

# Visualize
sc.pl.umap(adata, color=['leiden', 'CD3', 'CD8', 'CD20', 'panCK'])
```

### 空间可视化

在空间环境中可视化细胞：

```python
import matplotlib.pyplot as plt

# Spatial scatter plot
fig, ax = plt.subplots(figsize=(15, 15))

# Color by cell type
cell_types = adata.obs['cell_type'].unique()
colors = plt.cm.tab10(np.linspace(0, 1, len(cell_types)))

for i, cell_type in enumerate(cell_types):
    mask = adata.obs['cell_type'] == cell_type
    coords = adata.obsm['spatial'][mask]
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=[colors[i]],
        label=cell_type,
        s=5,
        alpha=0.7
    )

ax.legend(markerscale=2)
ax.set_xlabel('X (pixels)')
ax.set_ylabel('Y (pixels)')
ax.set_title('Spatial Cell Type Distribution')
ax.axis('equal')
plt.tight_layout()
plt.show()
```

### 空间邻域分析

分析细胞邻域和相互作用：

```python
import squidpy as sq

# Calculate spatial neighborhood enrichment
sq.gr.spatial_neighbors(adata, coord_type='generic', spatial_key='spatial')

# Neighborhood enrichment test
sq.gr.nhood_enrichment(adata, cluster_key='cell_type')

# Visualize interaction matrix
sq.pl.nhood_enrichment(adata, cluster_key='cell_type')

# Co-occurrence score
sq.gr.co_occurrence(adata, cluster_key='cell_type')
sq.pl.co_occurrence(
    adata,
    cluster_key='cell_type',
    clusters=['CD8 T cell', 'Tumor'],
    figsize=(8, 8)
)
```

### 空间自相关

测试标记的空间聚类：

```python
# Moran's I spatial autocorrelation
sq.gr.spatial_autocorr(
    adata,
    mode='moran',
    genes=['CD3', 'CD8', 'PD1', 'PDL1', 'panCK']
)

# Visualize
results = adata.uns['moranI']
print(results.head())
```

## MERFISH 工作流程

### 加载 MERFISH 数据

```python
from pathml.core import MERFISHSlide

# Load MERFISH dataset
merfish_slide = MERFISHSlide(
    path='path/to/merfish_data',
    fov_size=2048,  # Field of view size
    microns_per_pixel=0.108
)
```

### MERFISH 处理

```python
from pathml.preprocessing import Pipeline, DecodeMERFISH, SegmentMIF

merfish_pipeline = Pipeline([
    # 1. Decode barcodes to genes
    DecodeMERFISH(
        codebook='path/to/codebook.csv',
        error_correction=True,
        distance_threshold=0.5
    ),

    # 2. Cell segmentation
    SegmentMIF(
        nuclear_channel='DAPI',
        cytoplasm_channel='polyT',  # poly(T) stain for cell boundaries
        model='mesmer'
    ),

    # 3. Assign transcripts to cells
    AssignTranscripts(
        segmentation_mask_name='cell_segmentation',
        transcript_coords='decoded_spots'
    )
])

merfish_pipeline.run(merfish_slide)

# Output: AnnData with gene counts per cell
gene_expression = merfish_slide.cell_data
```

## 质量控制

### 分割质量

```python
from pathml.utils import assess_segmentation_quality

# Check segmentation quality metrics
qc_metrics = assess_segmentation_quality(
    segmentation_mask,
    image,
    metrics=['cell_count', 'mean_cell_size', 'size_distribution']
)

print(f"Total cells: {qc_metrics['cell_count']}")
print(f"Mean cell size: {qc_metrics['mean_cell_size']:.1f} pixels")

# Visualize
import matplotlib.pyplot as plt
plt.hist(qc_metrics['cell_sizes'], bins=50)
plt.xlabel('Cell Size (pixels)')
plt.ylabel('Frequency')
plt.title('Cell Size Distribution')
plt.show()
```

### 标记表达 QC

```python
import scanpy as sc

# Load AnnData
adata = ad.read_h5ad('codex_dataset.h5ad')

# Calculate QC metrics
adata.obs['total_intensity'] = adata.X.sum(axis=1)
adata.obs['n_markers_detected'] = (adata.X > 0).sum(axis=1)

# Filter low-quality cells
adata = adata[adata.obs['total_intensity'] > 100, :]
adata = adata[adata.obs['n_markers_detected'] >= 3, :]

# Visualize
sc.pl.violin(adata, ['total_intensity', 'n_markers_detected'], multi_panel=True)
```

## 批处理

高效处理大型多参数数据集：

```python
from pathml.core import SlideDataset
from pathml.preprocessing import Pipeline
from dask.distributed import Client
import glob

# Start Dask cluster
client = Client(n_workers=16, threads_per_worker=2, memory_limit='8GB')

# Find all CODEX slides
slide_dirs = glob.glob('data/codex_slides/*/')

# Create dataset
codex_slides = [CODEXSlide(d, stain='IF') for d in slide_dirs]
dataset = SlideDataset(codex_slides)

# Run pipeline in parallel
dataset.run(
    codex_pipeline,
    distributed=True,
    client=client,
    scheduler='distributed'
)

# Save processed data
for i, slide in enumerate(dataset):
    slide.cell_data.write(f'processed/slide_{i}.h5ad')

client.close()
```

## 与其他工具集成

### 导出到空间分析工具

```python
# Export to Giotto
def export_to_giotto(adata, output_dir):
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Expression matrix
    pd.DataFrame(
        adata.X.T,
        index=adata.var_names,
        columns=adata.obs_names
    ).to_csv(f'{output_dir}/expression.csv')

    # Cell coordinates
    pd.DataFrame(
        adata.obsm['spatial'],
        columns=['x', 'y'],
        index=adata.obs_names
    ).to_csv(f'{output_dir}/spatial_locs.csv')

# Export to Seurat
def export_to_seurat(adata, output_file):
    adata.write_h5ad(output_file)
    # Read in R with: library(Seurat); ReadH5AD(output_file)
```

## 最佳实践

1. **分段的通道选择：**
   - 使用最亮、最一致的核标记（通常是 DAPI）
   - 根据组织类型选择膜/细胞质标记
   - 测试多个选项以优化细分

2. **背景扣除：**
- 在定量前应用以减少自发荧光
   - 使用空白/控制图像来模拟背景

3. **质量控制：**
   - 可视化样本区域的分割
   - 检查像元大小分布中的异常值
   - 验证标记表达范围

4. **细胞类型注释：**
   - 从规范标记开始（CD3、CD20、panCK）
   - 使用多个标记进行稳健分类
   - 考虑无监督聚类来发现群体

5. **空间分析：**
   - 考虑组织结构（上皮、间质等）
   - 解释相互作用时考虑局部密度
   - 使用排列检验来确定统计显着性

6. **批量效果：**
   - 在 AnnData.obs 中包含批次信息
   - 如果组合多个实验，则应用批量校正
   - 使用按批次着色的 UMAP 可视化批次效果

## 常见问题及解决方案

**问题：分割质量差**
- 验证核和细胞质通道是否正确指定
- 调整 image_resolution 参数以匹配实际分辨率
- 尝试不同的细胞质标记
- 手动调整最小/最大单元尺寸参数

**问题：标记强度低**
- 检查背景扣除伪影
- 验证频道名称与实际频道匹配
- 检查原始图像是否存在技术问题（对焦、曝光）

**问题：细胞类型注释与预期不符**
- 调整标记阈值（太高/太低）
- 可视化标记分布以设置数据驱动的阈值
- 检查抗体特异性问题

**问题：空间分析显示没有显着的相互作用**
- 增加邻里半径
- 检查每种类型的细胞数量是否充足
- 验证空间坐标是否正确缩放

## 其他资源

- **PathML 多参数 API：** https://pathml.readthedocs.io/en/latest/api_multiparametric_reference.html
- **法典：** https://www.akoyabio.com/codex/
- **Vectra：** https://www.akoyabio.com/vectra/
- **DeepCell Mesmer：** https://www.deepcell.org/
- **Scanpy:** https://scanpy.readthedocs.io/（单细胞分析）
- **Squidpy:** https://squidpy.readthedocs.io/（空间组学分析）