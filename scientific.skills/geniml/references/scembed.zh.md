<!-- 此文件由机器翻译自 scembed.md -->

# scEmbed：单细胞嵌入生成

## 概述

scEmbed 在单细胞 ATAC-seq 数据集上训练 Region2Vec 模型，以生成用于聚类和分析的细胞嵌入。它提供了一个无监督机器学习框架，用于表示和分析低维空间中的 scATAC-seq 数据。

## 何时使用

在处理以下内容时使用 scEmbed：
- 需要聚类的单细胞 ATAC-seq (scATAC-seq) 数据
- 细胞类型注释任务
- 单细胞染色质可及性的降维
- 与 scanpy 工作流程集成以进行下游分析

## 工作流程

### 第 1 步：数据准备

输入数据必须采用 AnnData 格式，其中 `.var` 属性包含峰值的 `chr`、`start` 和 `end` 值。

**从原始数据开始**（barcodes.txt、peaks.bed、matrix.mtx）：

```python
import scanpy as sc
import pandas as pd
import scipy.io
import anndata

# Load data
barcodes = pd.read_csv('barcodes.txt', header=None, names=['barcode'])
peaks = pd.read_csv('peaks.bed', sep='\t', header=None,
                    names=['chr', 'start', 'end'])
matrix = scipy.io.mmread('matrix.mtx').tocsr()

# Create AnnData
adata = anndata.AnnData(X=matrix.T, obs=barcodes, var=peaks)
adata.write('scatac_data.h5ad')
```

### 第 2 步：预标记化

使用 gtars 实用程序将基因组区域转换为标记。这将创建一个包含标记化单元的镶木地板文件，以加快训练速度：

<<<代码块_1>>>

**预标记化的好处：**
- 更快的训练迭代
- 减少内存需求
- 可重复使用标记化数据进行多次训练

### 步骤 3：模型训练

使用标记化数据训练 scEmbed 模型：

<<<代码块_2>>>

### 步骤 4：生成单元嵌入

使用经过训练的模型生成单元格的嵌入：

<<<代码块_3>>>

### 步骤 5：下游分析

与 scanpy 集成以进行聚类和可视化：

<<<代码块_4>>>

## 关键参数

### 训练参数

|参数|描述 |典型范围|
|------------|-------------|---------------|
| `embedding_dim` |细胞嵌入的尺寸 | 50 - 200 | 50 - 200
| `window_size` |训练上下文窗口 | 3 - 10 | 3 - 10
| `negative_samples` |负样本数量| 5 - 20 | 5 - 20
| `epochs` |训练纪元 | 50 - 200 | 50 - 200
| `batch_size` |训练批量大小 | 128 - 512 | 128 - 512
| `learning_rate` |初始学习率| 0.01 - 0.05 | 0.01 - 0.05

### 标记化参数

- **Universe 文件**：定义基因组词汇的参考 BED 文件
- **重叠阈值**：峰宇宙匹配的最小重叠（通常为 1e-9）

## 预训练模型

Hugging Face 上提供了预训练的 scEmbed 模型，用于通用参考数据集。使用以下方式加载它们：

<<<代码块_5>>>

## 最佳实践

- **数据质量**：使用过滤后的峰值条形码矩阵，而不是原始计数
- **预分词**：始终进行预分词以提高训练效率
- **参数调整**：根据数据集大小调整`embedding_dim`和训练时期
- **验证**：使用已知的细胞类型标记来验证聚类质量
- **集成**：与scanpy结合进行全面的单细胞分析
- **模型共享**：将训练好的模型导出到 Hugging Face 以实现可重复性

## 示例数据集

10x Genomics PBMC 10k 数据集（10,000 个外周血单核细胞）作为标准基准：
- 含有多种免疫细胞类型
- 充分表征的细胞群
- 可从 10x Genomics 网站获取

## 单元格类型注释

聚类后，使用 k 最近邻 (KNN) 和参考数据集注释细胞类型：

<<<代码块_6>>>

## 输出

scEmbed 产生：
- 低维单元嵌入（存储在`adata.obsm`中）
- 经过训练的模型文件以供重复使用
- scanpy下游分析的兼容格式
- 可选择导出至 Hugging Face 以进行共享