<!-- 此文件由机器翻译自 plotting_guide.md -->

# Scanpy 绘图指南

使用 scanpy 创建出版质量可视化的综合指南。

## 一般绘图原则

所有 scanpy 绘图函数都遵循一致的模式：
- `sc.pl.*` 中的函数镜像 `sc.tl.*` 中的分析函数
- 大多数接受基因名称或元数据列的 `color` 参数
- 结果通过`save`参数保存
- 一次调用可以生成多个绘图

## 基本质量控制图

### 可视化 QC 指标

```python
# Violin plots for QC metrics
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save='_qc_violin.pdf')

# Scatter plots to identify outliers
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_qc_mt.pdf')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_qc_genes.pdf')

# Highest expressing genes
sc.pl.highest_expr_genes(adata, n_top=20, save='_highest_expr.pdf')
```

### 过滤后 QC

<<<代码块_1>>>

## 降维可视化

### 主成分分析图

<<<代码块_2>>>

### UMAP 图

<<<代码块_3>>>

### t-SNE 图

<<<代码块_4>>>

## 聚类可视化

### 基本聚类图

<<<代码块_5>>>

### 集群比较

<<<代码块_6>>>

## 标记基因可视化

### 排名标记基因

```python
# Overview of top markers per cluster
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,
                        save='_marker_overview.pdf')

# Heatmap of top markers
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby='leiden',
                                 show_gene_labels=True,
                                 save='_marker_heatmap.pdf')

# Dot plot of markers
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5,
                                 save='_marker_dotplot.pdf')

# Stacked violin plots
sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5,
                                        save='_marker_violin.pdf')

# Matrix plot
sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5,
                                    save='_marker_matrix.pdf')
```

### 特定基因表达

```python
# Violin plots for specific genes
marker_genes = ['CD3D', 'CD14', 'MS4A1', 'NKG7', 'FCGR3A']
sc.pl.violin(adata, keys=marker_genes, groupby='leiden',
             save='_markers_violin.pdf')

# Dot plot for curated markers
sc.pl.dotplot(adata, var_names=marker_genes, groupby='leiden',
              save='_markers_dotplot.pdf')

# Heatmap for specific genes
sc.pl.heatmap(adata, var_names=marker_genes, groupby='leiden',
              swap_axes=True, save='_markers_heatmap.pdf')

# Stacked violin for gene sets
sc.pl.stacked_violin(adata, var_names=marker_genes, groupby='leiden',
                     save='_markers_stacked.pdf')
```

### 嵌入的基因表达

```python
# Multiple genes on UMAP
genes = ['CD3D', 'CD14', 'MS4A1', 'NKG7']
sc.pl.umap(adata, color=genes, cmap='viridis',
           save='_umap_markers.pdf')

# Gene expression with custom colormap
sc.pl.umap(adata, color='CD3D', cmap='Reds',
           vmin=0, vmax=3, save='_umap_cd3d.pdf')
```

## 轨迹和伪时间可视化

### PAGA 情节

```python
# PAGA graph
sc.pl.paga(adata, color='leiden', save='_paga.pdf')

# PAGA with gene expression
sc.pl.paga(adata, color=['leiden', 'dpt_pseudotime'],
           save='_paga_pseudotime.pdf')

# PAGA overlaid on UMAP
sc.pl.umap(adata, color='leiden', save='_umap_with_paga.pdf',
           edges=True, edges_color='gray')
```

### 伪时间图

```python
# DPT pseudotime on UMAP
sc.pl.umap(adata, color='dpt_pseudotime', save='_umap_dpt.pdf')

# Gene expression along pseudotime
sc.pl.dpt_timeseries(adata, save='_dpt_timeseries.pdf')

# Heatmap ordered by pseudotime
sc.pl.heatmap(adata, var_names=genes, groupby='leiden',
              use_raw=False, show_gene_labels=True,
              save='_pseudotime_heatmap.pdf')
```

## 高级可视化

### 轨迹图（基因表达趋势）

```python
# Show gene expression across cell types
sc.pl.tracksplot(adata, var_names=marker_genes, groupby='leiden',
                 save='_tracks.pdf')
```

### 相关矩阵

```python
# Correlation between clusters
sc.pl.correlation_matrix(adata, groupby='leiden',
                         save='_correlation.pdf')
```

### 嵌入密度

```python
# Cell density on UMAP
sc.tl.embedding_density(adata, basis='umap', groupby='cell_type')
sc.pl.embedding_density(adata, basis='umap', key='umap_density_cell_type',
                        save='_density.pdf')
```

## 多面板人物

### 创建面板图

```python
import matplotlib.pyplot as plt

# Create multi-panel figure
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

# Plot on specific axes
sc.pl.umap(adata, color='leiden', ax=axes[0, 0], show=False)
sc.pl.umap(adata, color='CD3D', ax=axes[0, 1], show=False)
sc.pl.umap(adata, color='CD14', ax=axes[1, 0], show=False)
sc.pl.umap(adata, color='MS4A1', ax=axes[1, 1], show=False)

plt.tight_layout()
plt.savefig('figures/multi_panel.pdf')
plt.show()
```

## 出版物质量定制

### 高品质设置

```python
# Set publication-quality defaults
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(5, 5),
                               facecolor='white')

# Vector graphics output
sc.settings.figdir = './figures/'
sc.settings.file_format_figs = 'pdf'  # or 'svg'
```

### 自定义调色板

```python
# Use custom colors
custom_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
sc.pl.umap(adata, color='leiden', palette=custom_colors,
           save='_custom_colors.pdf')

# Continuous color maps
sc.pl.umap(adata, color='CD3D', cmap='viridis', save='_viridis.pdf')
sc.pl.umap(adata, color='CD3D', cmap='RdBu_r', save='_rdbu.pdf')
```

### 删除轴和框架

```python
# Clean plot without axes
sc.pl.umap(adata, color='leiden', frameon=False,
           save='_clean.pdf')

# No legend
sc.pl.umap(adata, color='leiden', legend_loc=None,
           save='_no_legend.pdf')
```

## 导出绘图

### 保存单独的图

```python
# Automatic saving with save parameter
sc.pl.umap(adata, color='leiden', save='_leiden.pdf')
# Saves to: sc.settings.figdir + 'umap_leiden.pdf'

# Manual saving
import matplotlib.pyplot as plt
fig = sc.pl.umap(adata, color='leiden', show=False, return_fig=True)
fig.savefig('figures/my_umap.pdf', dpi=300, bbox_inches='tight')
```

### 批量导出

```python
# Save multiple versions
for gene in ['CD3D', 'CD14', 'MS4A1']:
    sc.pl.umap(adata, color=gene, save=f'_{gene}.pdf')
```

## 常用自定义参数

### 布局参数
- `figsize`：图形尺寸（宽度、高度）
- `frameon`：在图周围显示框架
- `title`：绘图标题
- `legend_loc`：“右边距”、“数据上”、“最佳”或“无”
- `legend_fontsize`：图例的字体大小
- `size`：点大小

### 颜色参数
- `color`：用于着色的变量
- `palette`：调色板（例如，“Set1”、“viridis”）
- `cmap`：连续变量的颜色图
- `vmin`、`vmax`：色阶限制
- `use_raw`：使用基因表达的原始计数

### 保存参数
- `save`：保存的文件名后缀
- `show`：是否显示绘图
- `dpi`：光栅格式的分辨率

## 出版数据的提示

1. **使用矢量格式**：PDF或SVG用于可缩放图形
2. **高DPI**：光栅图像设置dpi=300或更高
3. **一致的样式**：在各个图形中使用相同的调色板
4. **清晰的标签**：确保基因名称和细胞类型可读
5. **白色背景**：使用 `facecolor='white'` 进行发布
6. **消除混乱**：设置 `frameon=False` 以获得更清晰的外观
7. **图例放置**：对于紧凑的图形使用“数据”
8. **色盲友好**：考虑像“colorblind”或“Set2”这样的调色板