<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 扫描
描述：“单细胞 RNA-seq 分析。加载 .h5ad/10X 数据、QC、标准化、PCA/UMAP/t-SNE、Leiden 聚类、标记基因、细胞类型注释、轨迹，用于 scRNA-seq 分析。”
---

# Scanpy：单细胞分析

## 概述

Scanpy 是一个可扩展的 Python 工具包，用于分析单细胞 RNA-seq 数据，基于 AnnData 构建。将此技能应用于完整的单细胞工作流程，包括质量控制、标准化、降维、聚类、标记基因识别、可视化和轨迹分析。

## 何时使用此技能

该技能应该在以下情况下使用：
- 分析单细胞 RNA-seq 数据（.h5ad、10X、CSV 格式）
- 对 scRNA-seq 数据集进行质量控制
- 创建 UMAP、t-SNE 或 PCA 可视化
- 识别细胞簇并寻找标记基因
- 根据基因表达注释细胞类型
- 进行轨迹推断或伪时间分析
- 生成出版质量的单细胞图

## 快速入门

### 基本导入和设置

```python
import scanpy as sc
import pandas as pd
import numpy as np

# Configure settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')
sc.settings.figdir = './figures/'
```

### 加载数据

<<<代码块_1>>>

### 了解 Ann 数据结构

AnnData对象是scanpy中的核心数据结构：

<<<代码块_2>>>

## 标准分析工作流程

### 1. 质量控制

识别并过滤低质量的细胞和基因：

<<<代码块_3>>>

**使用QC脚本进行自动分析：**
<<<代码块_4>>>

### 2. 标准化和预处理

<<<代码块_5>>>

### 3.降维

<<<代码块_6>>>

### 4. 聚类

```python
# Leiden clustering (recommended)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color='leiden', legend_loc='on data')

# Try multiple resolutions to find optimal granularity
for res in [0.3, 0.5, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
```

### 5. 标记基因鉴定

```python
# Find marker genes for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Visualize results
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)

# Get results as DataFrame
markers = sc.get.rank_genes_groups_df(adata, group='0')
```

### 6. 单元格类型注释

```python
# Define marker genes for known cell types
marker_genes = ['CD3D', 'CD14', 'MS4A1', 'NKG7', 'FCGR3A']

# Visualize markers
sc.pl.umap(adata, color=marker_genes, use_raw=True)
sc.pl.dotplot(adata, var_names=marker_genes, groupby='leiden')

# Manual annotation
cluster_to_celltype = {
    '0': 'CD4 T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    '3': 'CD8 T cells',
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype)

# Visualize annotated types
sc.pl.umap(adata, color='cell_type', legend_loc='on data')
```

### 7. 保存结果

```python
# Save processed data
adata.write('results/processed_data.h5ad')

# Export metadata
adata.obs.to_csv('results/cell_metadata.csv')
adata.var.to_csv('results/gene_metadata.csv')
```

## 常见任务

### 创建出版质量的绘图

```python
# Set high-quality defaults
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(5, 5))
sc.settings.file_format_figs = 'pdf'

# UMAP with custom styling
sc.pl.umap(adata, color='cell_type',
           palette='Set2',
           legend_loc='on data',
           legend_fontsize=12,
           legend_fontoutline=2,
           frameon=False,
           save='_publication.pdf')

# Heatmap of marker genes
sc.pl.heatmap(adata, var_names=genes, groupby='cell_type',
              swap_axes=True, show_gene_labels=True,
              save='_markers.pdf')

# Dot plot
sc.pl.dotplot(adata, var_names=genes, groupby='cell_type',
              save='_dotplot.pdf')
```

有关全面的可视化示例，请参阅`references/plotting_guide.md`。

### 轨迹推断

```python
# PAGA (Partition-based graph abstraction)
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color='leiden')

# Diffusion pseudotime
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == '0')[0]
sc.tl.dpt(adata)
sc.pl.umap(adata, color='dpt_pseudotime')
```

### 条件之间的差异表达式

```python
# Compare treated vs control within cell types
adata_subset = adata[adata.obs['cell_type'] == 'T cells']
sc.tl.rank_genes_groups(adata_subset, groupby='condition',
                         groups=['treated'], reference='control')
sc.pl.rank_genes_groups(adata_subset, groups=['treated'])
```

### 基因集评分

```python
# Score cells for gene set expression
gene_set = ['CD3D', 'CD3E', 'CD3G']
sc.tl.score_genes(adata, gene_set, score_name='T_cell_score')
sc.pl.umap(adata, color='T_cell_score')
```

### 批量修正

```python
# ComBat batch correction
sc.pp.combat(adata, key='batch')

# Alternative: use Harmony or scVI (separate packages)
```

## 需要调整的关键参数

### 质量控制
- `min_genes`：每个细胞的最小基因数（通常为 200-500）
- `min_cells`：每个基因的最小细胞数（通常为 3-10）
- `pct_counts_mt`：线粒体阈值（通常为 5-20%）

### 标准化
- `target_sum`：每个单元格的目标计数（默认 1e4）

### 特征选择
- `n_top_genes`：HVG 数量（通常为 2000-3000）
- `min_mean`、`max_mean`、`min_disp`：HVG 选择参数

### 降维
- `n_pcs`：主成分数量（检查方差比图）
- `n_neighbors`：邻居数量（通常为 10-30）

### 聚类
- `resolution`：聚类粒度（0.4-1.2，更高=更多聚类）

## 常见陷阱和最佳实践

1. **始终保存原始计数**：在过滤基因之前`adata.raw = adata`
2. **仔细检查 QC 图**：根据数据集质量调整阈值
3. **使用莱顿而不是鲁汶**：更高效、效果更好
4. **尝试多种聚类分辨率**：找到最佳粒度
5. **验证细胞类型注释**：使用多个标记基因
6. **使用 `use_raw=True` 绘制基因表达图**：显示原始计数
7. **检查PCA方差比**：确定最佳PC数量
8. **保存中间结果**：长工作流程可能会中途失败

## 捆绑资源

### 脚本/qc_analysis.py
自动质量控制脚本，用于计算指标、生成绘图和过滤数据：

```bash
python scripts/qc_analysis.py input.h5ad --output filtered.h5ad \
    --mt-threshold 5 --min-genes 200 --min-cells 3
```

### 参考文献/standard_workflow.md
完整的分步工作流程，包含以下详细说明和代码示例：
- 数据加载和设置
- 可视化质量控制
- 标准化和缩放
- 特征选择
- 降维（PCA、UMAP、t-SNE）
- 聚类（莱顿、鲁汶）
- 标记基因鉴定
- 细胞类型注释
- 轨迹推断
- 差异表达

从头开始执行完整分析时请阅读此参考资料。

### 参考资料/api_reference.md
按模块组织的 scanpy 函数快速参考指南：
- 读取/写入数据（`sc.read_*`、`adata.write_*`）
- 预处理（`sc.pp.*`）
- 工具（`sc.tl.*`）
- 绘图（`sc.pl.*`）
- Ann数据结构和操作
- 设置和实用程序

使用它可以快速查找函数签名和通用参数。

### 参考文献/plotting_guide.md
全面的可视化指南包括：
- 质量控制图
- 降维可视化
- 聚类可视化
- 标记基因图（热图、点图、小提琴图）
- 轨迹和伪时间图
- 出版物质量的定制
- 多面板人物
- 调色板和造型

创建可供出版的图表时请参考此内容。

### 资产/analysis_template.py
完整的分析模板提供从数据加载到细胞类型注释的完整工作流程。复制并自定义此模板以进行新分析：

```bash
cp assets/analysis_template.py my_analysis.py
# Edit parameters and run
python my_analysis.py
```

该模板包括所有标准步骤以及可配置参数和有用的注释。

## 其他资源

- **官方 scanpy 文档**：https://scanpy.readthedocs.io/
- **Scanpy 教程**：https://scanpy-tutorials.readthedocs.io/
- **scverse 生态系统**：https://scverse.org/（相关工具：squidpy、scvi-tools、cellrank）
- **最佳实践**：Luecken & Theis (2019)“单细胞 RNA 测序的当前最佳实践”

## 有效分析的技巧

1. **从模板开始**：使用 `assets/analysis_template.py` 作为起点
2. **首先运行QC脚本**：使用`scripts/qc_analysis.py`进行初始过滤
3. **根据需要查阅参考资料**：将工作流程和 API 参考加载到上下文中
4. **迭代聚类**：尝试多种分辨率和可视化方法
5. **生物学验证**：检查标记基因是否与预期细胞类型相匹配
6. **文件参数**：记录QC阈值和分析设置
7. **保存检查点**：在关键步骤写入中间结果