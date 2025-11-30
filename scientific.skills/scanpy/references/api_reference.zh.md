<!-- 此文件由机器翻译自 api_reference.md -->

# Scanpy API 快速参考

按模块组织的常用 scanpy 函数的快速参考。

## 进口公约

```python
import scanpy as sc
```

## 读取和写入数据 (sc.read_*)

### 读取函数

<<<代码块_1>>>

### 编写函数

<<<代码块_2>>>

## 预处理 (sc.pp.*)

### 质量控制

<<<代码块_3>>>

### 标准化和转换

<<<代码块_4>>>

### 特征选择

<<<代码块_5>>>

### 缩放和回归

<<<代码块_6>>>

### 降维（预处理）

```python
sc.pp.pca(adata, n_comps=50)                     # Principal component analysis
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40) # Compute neighborhood graph
```

### 批量修正

```python
sc.pp.combat(adata, key='batch')                 # ComBat batch correction
```

## 工具 (sc.tl.*)

### 降维

```python
sc.tl.pca(adata, svd_solver='arpack')            # PCA
sc.tl.umap(adata)                                 # UMAP embedding
sc.tl.tsne(adata)                                 # t-SNE embedding
sc.tl.diffmap(adata)                              # Diffusion map
sc.tl.draw_graph(adata, layout='fa')             # Force-directed graph
```

### 聚类

```python
sc.tl.leiden(adata, resolution=0.5)              # Leiden clustering (recommended)
sc.tl.louvain(adata, resolution=0.5)             # Louvain clustering
sc.tl.kmeans(adata, n_clusters=10)               # K-means clustering
```

### 标记基因和差异表达

```python
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
sc.tl.rank_genes_groups(adata, groupby='leiden', method='logreg')

# Get results as dataframe
sc.get.rank_genes_groups_df(adata, group='0')
```

### 轨迹推断

```python
sc.tl.paga(adata, groups='leiden')               # PAGA trajectory
sc.tl.dpt(adata)                                  # Diffusion pseudotime
```

### 基因评分

```python
sc.tl.score_genes(adata, gene_list, score_name='score')
sc.tl.score_genes_cell_cycle(adata, s_genes, g2m_genes)
```

### 嵌入和投影

```python
sc.tl.ingest(adata, adata_ref)                   # Map to reference
sc.tl.embedding_density(adata, basis='umap', groupby='leiden')
```

## 绘图 (sc.pl.*)

### 基本嵌入

```python
sc.pl.umap(adata, color='leiden')                # UMAP plot
sc.pl.tsne(adata, color='gene_name')             # t-SNE plot
sc.pl.pca(adata, color='leiden')                 # PCA plot
sc.pl.diffmap(adata, color='leiden')             # Diffusion map plot
```

### 热图和点图

```python
sc.pl.heatmap(adata, var_names=genes, groupby='leiden')
sc.pl.dotplot(adata, var_names=genes, groupby='leiden')
sc.pl.matrixplot(adata, var_names=genes, groupby='leiden')
sc.pl.stacked_violin(adata, var_names=genes, groupby='leiden')
```

### 小提琴图和散点图

```python
sc.pl.violin(adata, keys=['gene1', 'gene2'], groupby='leiden')
sc.pl.scatter(adata, x='gene1', y='gene2', color='leiden')
```

### 标记基因可视化

```python
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
sc.pl.rank_genes_groups_violin(adata, groups='0')
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)
```

### 轨迹可视化

```python
sc.pl.paga(adata, color='leiden')                # PAGA graph
sc.pl.dpt_timeseries(adata)                      # DPT timeseries
```

### QC 图

```python
sc.pl.highest_expr_genes(adata, n_top=20)
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'])
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
```

### 高级情节

```python
sc.pl.dendrogram(adata, groupby='leiden')
sc.pl.correlation_matrix(adata, groupby='leiden')
sc.pl.tracksplot(adata, var_names=genes, groupby='leiden')
```

## 常用参数

### 颜色参数
- `color`：用于着色的变量（基因名称、obs 列）
- `use_raw`：使用adata的`.raw`属性
- `palette`：要使用的调色板
- `vmin`、`vmax`：色阶限制

### 布局参数
- `basis`：嵌入基础（'umap'、'tsne'、'pca' 等）
- `legend_loc`：图例位置（“数据上”、“右边距”等）
- `size`：点大小
- `alpha`：点透明度

### 保存参数
- `save`：保存绘图的文件名
- `show`：是否显示绘图

## Ann数据结构

```python
adata.X                    # Expression matrix (cells × genes)
adata.obs                  # Cell annotations (DataFrame)
adata.var                  # Gene annotations (DataFrame)
adata.uns                  # Unstructured annotations (dict)
adata.obsm                 # Multi-dimensional cell annotations (e.g., PCA, UMAP)
adata.varm                 # Multi-dimensional gene annotations
adata.layers               # Additional data layers
adata.raw                  # Raw data backup

# Access
adata.obs_names            # Cell barcodes
adata.var_names            # Gene names
adata.shape                # (n_cells, n_genes)

# Slicing
adata[cell_indices, gene_indices]
adata[:, adata.var_names.isin(gene_list)]
adata[adata.obs['leiden'] == '0', :]
```

## 设置

```python
sc.settings.verbosity = 3              # 0=error, 1=warning, 2=info, 3=hint
sc.settings.set_figure_params(dpi=80, facecolor='white')
sc.settings.autoshow = False           # Don't show plots automatically
sc.settings.autosave = True            # Autosave figures
sc.settings.figdir = './figures/'      # Figure directory
sc.settings.cachedir = './cache/'      # Cache directory
sc.settings.n_jobs = 8                 # Number of parallel jobs
```

## 有用的实用程序

```python
sc.logging.print_versions()            # Print version information
sc.logging.print_memory_usage()        # Print memory usage
adata.copy()                           # Create a copy of AnnData object
adata.concatenate([adata1, adata2])    # Concatenate AnnData objects
```