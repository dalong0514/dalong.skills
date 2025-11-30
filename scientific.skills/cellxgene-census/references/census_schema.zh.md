<!-- 此文件由机器翻译自 census_schema.md -->

# CZ CELLxGENE 人口普查数据架构参考

## 概述

CZ CELLxGENE Census 是基于 TileDB-SOMA 框架构建的单细胞数据的版本集合。此参考记录了数据结构、可用元数据字段和查询语法。

## 高层结构

人口普查以 `SOMACollection` 的形式组织，具有两个主要组成部分：

### 1.人口普查信息
摘要信息包括：
- **摘要**：构建日期、细胞计数、数据集统计信息
- **数据集**：来自 CELLxGENE Discover 的所有数据集以及元数据
- **summary_cell_counts**：按元数据类别分层的细胞计数

### 2. 人口普查数据
生物体特定的 `SOMAExperiment` 对象：
- **“homo_sapiens”**：人类单细胞数据
- **“mus_musculus”**：小鼠单细胞数据

## 每个生物体的数据结构

每个生物体实验包含：

### obs（细胞元数据）
单元格级注释存储为 `SOMADataFrame`。访问方式：
```python
census["census_data"]["homo_sapiens"].obs
```

### ms["RNA"]（测量）
RNA测量数据包括：
- **X**：带层的数据矩阵：
  - `raw`：原始计数数据
  - `normalized`：（如果可用）标准化计数
- **var**：基因元数据
- **feature_dataset_presence_matrix**：稀疏布尔数组显示每个数据集中测量了哪些基因

## 单元格元数据字段 (obs)

### 必填/核心字段

**身份和数据集：**
- `soma_joinid`：连接的唯一整数标识符
- `dataset_id`：源数据集标识符
- `is_primary_data`：布尔标志（True = 唯一单元格，False = 跨数据集重复）

**细胞类型：**
- `cell_type`：人类可读的单元格类型名称
- `cell_type_ontology_term_id`：标准化本体术语（例如“CL:0000236”）

**纸巾：**
- `tissue`：特定组织名称
- `tissue_general`：更广泛的组织类别（对于分组有用）
- `tissue_ontology_term_id`：标准化本体术语

**化验：**
- `assay`：使用的测序技术
- `assay_ontology_term_id`：标准化本体术语

**疾病：**
- `disease`：疾病状态或状况
- `disease_ontology_term_id`：标准化本体术语

**捐助者：**
- `donor_id`：唯一的捐赠者标识符
- `sex`：生物性别（男、女、未知）
- `self_reported_ethnicity`：种族信息
- `development_stage`：生命阶段（成人、儿童、胚胎等）
- `development_stage_ontology_term_id`：标准化本体术语

**生物体：**
- `organism`：学名（Homo sapiens，Mus musculus）
- `organism_ontology_term_id`：标准化本体术语

**技术：**
- `suspension_type`：样品制备类型（细胞、细胞核、na）

## 基因元数据字段 (var)

访问方式：
<<<代码块_1>>>

**可用字段：**
- `soma_joinid`：连接的唯一整数标识符
- `feature_id`：Ensembl 基因 ID（例如“ENSG00000161798”）
- `feature_name`：基因符号（例如“FOXP2”）
- `feature_length`：碱基对中的基因长度

## 值过滤器语法

查询使用类似 Python 的表达式进行过滤。语法由 TileDB-SOMA 处理。

### 比较运算符
- `==`：等于
- `!=`：不等于
- `<`、`>`、`<=`、`>=`：数字比较
- `in`：成员资格测试（例如，`feature_id in ['ENSG00000161798', 'ENSG00000188229']`）

### 逻辑运算符
- `and`、`&`：逻辑与
- `or`、`|`：逻辑或

### 示例

**单一条件：**
<<<代码块_2>>>

**使用 AND 的多个条件：**
<<<代码块_3>>>

**使用 IN 表示多个值：**
<<<代码块_4>>>

**复杂情况：**
<<<代码块_5>>>

**过滤基因：**
<<<代码块_6>>>

## 数据包含标准

人口普查包括 CZ CELLxGENE Discover 会议的所有数据：

1. **物种**：人类 (*Homo sapiens*) 或小鼠 (*Mus musculus*)
2. **技术**：已批准的 RNA 测序技术
3. **计数类型**：仅原始计数（无处理/仅标准化数据）
4. **元数据**：标准化以下 CELLxGENE 模式
5. **空间和非空间数据**：包括传统和空间转录组学

## 重要数据特征

### 重复单元格
单元格可能出现在多个数据集中。在大多数分析中使用 `is_primary_data == True` 过滤唯一的单元格。

### 计数类型
人口普查包括：
- **分子计数**：来自基于 UMI 的方法
- **全基因测序读取计数**：来自非 UMI 方法
这些可能需要不同的标准化方法。

### 版本控制
人口普查版本有版本号（例如“2023-07-25”、“稳定”）。始终指定可重复分析的版本：
```python
census = cellxgene_census.open_soma(census_version="2023-07-25")
```

## 数据集存在矩阵

访问每个数据集中测量的基因：
```python
presence_matrix = census["census_data"]["homo_sapiens"].ms["RNA"]["feature_dataset_presence_matrix"]
```

这个稀疏布尔矩阵有助于理解：
- 跨数据集的基因覆盖率
- 特定基因分析应包含哪些数据集
- 与基因覆盖率相关的技术批次效应

## SOMA 对象类型

使用的核心 TileDB-SOMA 对象：
- **DataFrame**：表格数据（obs、var）
- **SparseNDArray**：稀疏矩阵（X层，存在矩阵）
- **DenseNDArray**：密集数组（不太常见）
- **Collection**：相关对象的容器
- **实验**：用于测量的顶级容器
- **SOMAScene**：空间转录组场景
- **obs_spatial_presence**：空间数据可用性