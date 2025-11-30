<!-- 此文件由机器翻译自 common_patterns.md -->

# 常见查询模式和最佳实践

## 查询模式类别

### 1. 探索性查询（仅限元数据）

在探索可用数据而不加载表达矩阵时使用。

**模式：获得组织中独特的细胞类型**
```python
import cellxgene_census

with cellxgene_census.open_soma() as census:
    cell_metadata = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter="tissue_general == 'brain' and is_primary_data == True",
        column_names=["cell_type"]
    )
    unique_cell_types = cell_metadata["cell_type"].unique()
    print(f"Found {len(unique_cell_types)} unique cell types")
```

**模式：按条件对单元格进行计数**
<<<代码块_1>>>

**模式：探索数据集信息**
<<<代码块_2>>>

### 2.中小型查询（AnnData）

当结果适合内存时（通常 < 100k 单元），请使用 `get_anndata()`。

**模式：组织特异性细胞类型查询**
<<<代码块_3>>>

**模式：具有多个基因的基因特异性查询**
<<<代码块_4>>>

**模式：多组织查询**
<<<代码块_5>>>

**模式：疾病特定查询**
<<<代码块_6>>>

### 3. 大型查询（核外处理）

对于超出可用 RAM 的查询，请使用 `axis_query()`。

**模式：迭代处理**
```python
import pyarrow as pa

# Create query
query = census["census_data"]["homo_sapiens"].axis_query(
    measurement_name="RNA",
    obs_query=soma.AxisQuery(
        value_filter="tissue_general == 'brain' and is_primary_data == True"
    ),
    var_query=soma.AxisQuery(
        value_filter="feature_name in ['FOXP2', 'TBR1', 'SATB2']"
    )
)

# Iterate through X matrix in chunks
iterator = query.X("raw").tables()
for batch in iterator:
    # Process batch (a pyarrow.Table)
    # batch has columns: soma_data, soma_dim_0, soma_dim_1
    process_batch(batch)
```

**模式：增量统计（均值/方差）**
```python
# Using Welford's online algorithm
n = 0
mean = 0
M2 = 0

iterator = query.X("raw").tables()
for batch in iterator:
    values = batch["soma_data"].to_numpy()
    for x in values:
        n += 1
        delta = x - mean
        mean += delta / n
        delta2 = x - mean
        M2 += delta * delta2

variance = M2 / (n - 1) if n > 1 else 0
```

### 4. PyTorch 集成（机器学习）

使用 `experiment_dataloader()` 来训练模型。

**模式：创建训练数据加载器**
```python
from cellxgene_census.experimental.ml import experiment_dataloader
import torch

with cellxgene_census.open_soma() as census:
    # Create dataloader
    dataloader = experiment_dataloader(
        census["census_data"]["homo_sapiens"],
        measurement_name="RNA",
        X_name="raw",
        obs_value_filter="tissue_general == 'liver' and is_primary_data == True",
        obs_column_names=["cell_type"],
        batch_size=128,
        shuffle=True,
    )

    # Training loop
    for epoch in range(num_epochs):
        for batch in dataloader:
            X = batch["X"]  # Gene expression
            labels = batch["obs"]["cell_type"]  # Cell type labels
            # Train model...
```

**模式：训练/测试分割**
```python
from cellxgene_census.experimental.ml import ExperimentDataset

# Create dataset from query
dataset = ExperimentDataset(
    experiment_axis_query,
    layer_name="raw",
    obs_column_names=["cell_type"],
    batch_size=128,
)

# Split data
train_dataset, test_dataset = dataset.random_split(
    split=[0.8, 0.2],
    seed=42
)

# Create loaders
train_loader = experiment_dataloader(train_dataset)
test_loader = experiment_dataloader(test_dataset)
```

### 5. 集成工作流程

**模式：Scanpy 集成**
```python
import scanpy as sc

# Load data
adata = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    obs_value_filter="cell_type == 'neuron' and is_primary_data == True",
)

# Standard scanpy workflow
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=["cell_type", "tissue_general"])
```

**模式：多数据集集成**
```python
# Query multiple datasets separately
datasets_to_integrate = ["dataset_id_1", "dataset_id_2", "dataset_id_3"]

adatas = []
for dataset_id in datasets_to_integrate:
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",
        obs_value_filter=f"dataset_id == '{dataset_id}' and is_primary_data == True",
    )
    adatas.append(adata)

# Integrate using scanorama, harmony, or other tools
import scanpy.external as sce
sce.pp.scanorama_integrate(adatas)
```

## 最佳实践

### 1. 始终过滤主要数据
除非专门分析重复项，否则始终包含 `is_primary_data == True`：
```python
obs_value_filter="cell_type == 'B cell' and is_primary_data == True"
```

### 2.指定人口普查版本
对于可重复的分析，请始终指定人口普查版本：
```python
census = cellxgene_census.open_soma(census_version="2023-07-25")
```

### 3.使用上下文管理器
始终使用上下文管理器来确保正确的清理：
```python
with cellxgene_census.open_soma() as census:
    # Your code here
```

### 4. 仅选择需要的列
通过仅选择所需的元数据列来最大限度地减少数据传输：
```python
obs_column_names=["cell_type", "tissue_general", "disease"]  # Not all columns
```

### 5. 检查基因查询的数据集是否存在
分析特定基因时，检查哪些数据集测量了它们：
```python
presence = cellxgene_census.get_presence_matrix(
    census,
    "homo_sapiens",
    var_value_filter="feature_name in ['CD4', 'CD8A']"
)
```

### 6.使用tissue_general进行更广泛的查询
`tissue_general` 提供比 `tissue` 更粗略的分组，对于跨组织分析非常有用：
```python
# Better for broad queries
obs_value_filter="tissue_general == 'immune system'"

# Use specific tissue when needed
obs_value_filter="tissue == 'peripheral blood mononuclear cell'"
```

### 7. 将元数据探索与表达式查询相结合
首先探索元数据以了解可用数据，然后查询表达式：
```python
# Step 1: Explore
metadata = cellxgene_census.get_obs(
    census, "homo_sapiens",
    value_filter="disease == 'COVID-19'",
    column_names=["cell_type", "tissue_general"]
)
print(metadata.value_counts())

# Step 2: Query based on findings
adata = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    obs_value_filter="disease == 'COVID-19' and cell_type == 'T cell' and is_primary_data == True",
)
```

### 8.大型查询的内存管理
对于大型查询，请在加载之前检查估计大小：
```python
# Get cell count first
metadata = cellxgene_census.get_obs(
    census, "homo_sapiens",
    value_filter="tissue_general == 'brain' and is_primary_data == True",
    column_names=["soma_joinid"]
)
n_cells = len(metadata)
print(f"Query will return {n_cells} cells")

# If too large, use out-of-core processing or further filtering
```

### 9. 利用本体术语来保持一致性
如果可能，请使用本体术语 ID 而不是自由文本：
```python
# More reliable than cell_type == 'B cell' across datasets
obs_value_filter="cell_type_ontology_term_id == 'CL:0000236'"
```

### 10.批处理模式
对于多种条件下的系统分析：
```python
tissues = ["lung", "liver", "kidney", "heart"]
results = {}

for tissue in tissues:
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="Homo sapiens",
        obs_value_filter=f"tissue_general == '{tissue}' and is_primary_data == True",
    )
    # Perform analysis
    results[tissue] = analyze(adata)
```

## 要避免的常见陷阱

1. **不过滤 is_primary_data**：导致计算重复单元格
2. **加载太多数据**：首先使用元数据查询来估计大小
3. **不使用上下文管理器**：可能导致资源泄漏
4. **版本控制不一致**：如果不指定版本，结果将无法重现
5. **过于宽泛的查询**：从有针对性的查询开始，根据需要扩展
6. **忽略数据集的存在**：某些基因未在所有数据集中测量
7. **错误的计数标准化**：注意 UMI 与读取计数的差异