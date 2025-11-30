<!-- 此文件由机器翻译自 standard_workflow.md -->

# 单细胞分析的标准 Scanpy 工作流程

本文档概述了使用 scanpy 分析单细胞 RNA-seq 数据的标准工作流程。

## 完整的分析流程

### 1. 数据加载和初始设置

```python
import scanpy as sc
import pandas as pd
import numpy as np

# Configure scanpy settings
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Load data (various formats)
adata = sc.read_10x_mtx('path/to/data/')  # For 10X data
# adata = sc.read_h5ad('path/to/data.h5ad')  # For h5ad format
# adata = sc.read_csv('path/to/data.csv')  # For CSV format
```

### 2. 质量控制 (QC)

<<<代码块_1>>>

### 3. 标准化

<<<代码块_2>>>

### 4.特征选择

<<<代码块_3>>>

### 5. 缩放和回归

<<<代码块_4>>>

### 6.降维

<<<代码块_5>>>

### 7. 聚类

<<<代码块_6>>>

### 8. 标记基因鉴定

```python
# Find marker genes for each cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Visualize top marker genes
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Get marker gene dataframe
marker_genes = sc.get.rank_genes_groups_df(adata, group='0')

# Visualize specific markers
sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
```

### 9. 单元格类型注释

```python
# Manual annotation based on marker genes
cluster_annotations = {
    '0': 'CD4 T cells',
    '1': 'CD14+ Monocytes',
    '2': 'B cells',
    '3': 'CD8 T cells',
    # ... add more annotations
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_annotations)

# Visualize annotated cell types
sc.pl.umap(adata, color='cell_type', legend_loc='on data')
```

### 10. 保存结果

```python
# Save the processed AnnData object
adata.write('results/processed_data.h5ad')

# Export results to CSV
adata.obs.to_csv('results/cell_metadata.csv')
adata.var.to_csv('results/gene_metadata.csv')
```

## 其他分析选项

### 轨迹推断

```python
# PAGA (Partition-based graph abstraction)
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color=['leiden'])

# Diffusion pseudotime (DPT)
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden'] == '0')[0]
sc.tl.dpt(adata)
sc.pl.umap(adata, color=['dpt_pseudotime'])
```

### 条件之间的差异表达式

```python
# Compare conditions within a cell type
sc.tl.rank_genes_groups(adata, groupby='condition', groups=['treated'],
                         reference='control', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['treated'])
```

### 基因集评分

```python
# Score cells for gene set expression
gene_set = ['CD3D', 'CD3E', 'CD3G']
sc.tl.score_genes(adata, gene_set, score_name='T_cell_score')
sc.pl.umap(adata, color='T_cell_score')
```

## 常用参数调整

- **QC 阈值**：`min_genes`、`min_cells`、`pct_counts_mt` - 取决于数据集质量
- **标准化目标**：通常为 1e4，但可以调整
- **HVG参数**：影响特征选择的严格性
- **PCA 成分**：检查方差比图以确定最佳数量
- **聚类分辨率**：值越高，聚类越多（通常为 0.4-1.2）
- **n_neighbors**：影响 UMAP 和聚类的粒度（通常为 10-30）

## 最佳实践

1. 在过滤之前始终可视化 QC 指标
2. 保存标准化前的原始计数 (`adata.raw = adata`)
3.使用Leiden代替Louvain进行聚类（效率更高）
4.尝试多种聚类分辨率以找到最佳粒度
5. 使用已知标记基因验证细胞类型注释
6.保存关键步骤的中间结果