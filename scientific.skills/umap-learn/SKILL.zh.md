<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：umap-学习
描述：“UMAP 降维。用于 2D/3D 可视化的快速非线性流形学习、聚类预处理 (HDBSCAN)、监督/参数化 UMAP，适用于高维数据。”
---

# UMAP-学习

## 概述

UMAP（统一流形逼近和投影）是一种用于可视化和一般非线性降维的降维技术。应用此技能来实现快速、可扩展的嵌入，从而保留局部和全局结构、监督学习和聚类预处理。

## 快速入门

### 安装

```bash
uv pip install umap-learn
```

### 基本用法

UMAP 遵循 scikit-learn 约定，可用作 t-SNE 或 PCA 的直接替代品。

<<<代码块_1>>>

**关键的预处理要求：** 在应用 UMAP 之前始终将特征标准化为可比较的尺度，以确保跨维度的权重相等。

### 典型工作流程

<<<代码块_2>>>

## 参数调优指南

UMAP 有四个控制嵌入行为的主要参数。了解这些对于有效使用至关重要。

### n_neighbors（默认值：15）

**目的：** 平衡嵌入中的局部结构与全局结构。

**工作原理：** 控制学习流形结构时 UMAP 检查的局部邻域的大小。

**按价值计算的影响：**
- **低值 (2-5)：** 强调精细的局部细节，但可能将数据分段为断开连接的组件
- **中等值 (15-20)：** 局部结构和全局关系的平衡视图（推荐起点）
- **高值 (50-200)：** 优先考虑广泛的拓扑结构，但牺牲细粒度的细节

**建议：** 从 15 开始，然后根据结果进行调整。增加以获得更多全局结构，减少以获得更多局部细节。

### min_dist（默认值：0.1）

**用途：**控制点在低维空间中聚集的紧密程度。

**工作原理：** 设置输出表示中允许点之间的最小距离。

**按价值计算的影响：**
- **低值 (0.0-0.1)：** 创建对聚类有用的聚集嵌入；揭示精细的拓扑细节
- **高值 (0.5-0.99)：** 防止紧密堆积；强调广泛的拓扑保存而不是局部结构

**建议：** 对于聚类应用程序使用 0.0，对于可视化使用 0.1-0.3，对于松散结构使用 0.5+。

### n_components（默认值：2）

**目的：**确定嵌入输出空间的维数。

**关键特征：** 与 t-SNE 不同，UMAP 在嵌入维度上可以很好地扩展，从而能够在可视化之外使用。

**常见用途：**
- **2-3 维度：** 可视化
- **5-10 维：** 聚类预处理（比 2D 更好地保留密度）
- **10-50 维度：** 下游 ML 模型的特征工程

**建议：** 使用 2 进行可视化，使用 5-10 进行聚类，使用更高的 ML 管道。

### 公制（默认值：'euclidean'）

**用途：** 指定如何计算输入数据点之间的距离。

**支持的指标：**
- **闵可夫斯基变体：** 欧几里得、曼哈顿、切比雪夫
- **空间度量：** 堪培拉、braycurtis、半正矢
- **相关性指标：**余弦、相关性（适用于文本/文档嵌入）
- **二进制数据指标：** Hamming、jaccard、dice、russellrao、kulsinski、rogerstanimoto、sokalmichener、sokalsneath、yule
- **自定义指标：** 通过 Numba 用户定义的距离函数

**建议：** 对于数值数据使用欧几里得，对于文本/文档向量使用余弦，对于二进制数据使用汉明。

### 参数调优示例

<<<代码块_3>>>

## 监督和半监督降维

UMAP 支持合并标签信息来指导嵌入过程，从而在保留内部结构的同时实现类分离。

### 监督UMAP

拟合时通过 `y` 参数传递目标标签：

<<<代码块_4>>>

**主要优点：**
- 实现清晰的班级划分
- 保留每个类的内部结构
- 维护类之间的全局关系

**何时使用：** 当您已标记数据并希望分离已知类同时保持有意义的点嵌入时。

### 半监督 UMAP

对于部分标签，请按照 scikit-learn 约定使用 `-1` 标记未标记的点：

<<<代码块_5>>>

**何时使用：** 当标签价格昂贵或您的数据多于可用标签时。

### 使用 UMAP 进行度量学习
在标记数据上训练监督嵌入，然后应用于新的未标记数据：

<<<代码块_6>>>

**何时使用：** 用于机器学习管道中的监督特征工程。

## 用于聚类的 UMAP

UMAP 可作为 HDBSCAN 等基于密度的聚类算法的有效预处理，克服维数灾难。

### 集群最佳实践

**关键原则：** 针对聚类和可视化配置 UMAP 的方式不同。

**推荐参数：**
- **n_neighbors：**增加到~30（默认15太本地化，可以创建人工细粒度集群）
- **min_dist：** 设置为 0.0（将点密集地打包在簇内以获得更清晰的边界）
- **n_components：** 使用 5-10 个维度（与 2D 相比，保持性能，同时提高密度保留）

### 聚类工作流程

```python
import umap
import hdbscan
from sklearn.preprocessing import StandardScaler

# 1. Preprocess data
scaled_data = StandardScaler().fit_transform(data)

# 2. UMAP with clustering-optimized parameters
reducer = umap.UMAP(
    n_neighbors=30,
    min_dist=0.0,
    n_components=10,  # Higher than 2 for better density preservation
    metric='euclidean',
    random_state=42
)
embedding = reducer.fit_transform(scaled_data)

# 3. Apply HDBSCAN clustering
clusterer = hdbscan.HDBSCAN(
    min_cluster_size=15,
    min_samples=5,
    metric='euclidean'
)
labels = clusterer.fit_predict(embedding)

# 4. Evaluate
from sklearn.metrics import adjusted_rand_score
score = adjusted_rand_score(true_labels, labels)
print(f"Adjusted Rand Score: {score:.3f}")
print(f"Number of clusters: {len(set(labels)) - (1 if -1 in labels else 0)}")
print(f"Noise points: {sum(labels == -1)}")
```

### 聚类后的可视化

```python
# Create 2D embedding for visualization (separate from clustering)
vis_reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
vis_embedding = vis_reducer.fit_transform(scaled_data)

# Plot with cluster labels
import matplotlib.pyplot as plt
plt.scatter(vis_embedding[:, 0], vis_embedding[:, 1], c=labels, cmap='Spectral', s=5)
plt.colorbar()
plt.title('UMAP Visualization with HDBSCAN Clusters')
plt.show()
```

**重要警告：** UMAP 不能完全保留密度，并且可以创建人为的簇划分。始终验证和探索生成的集群。

## 转换新数据

UMAP 通过其 `transform()` 方法启用新数据的预处理，允许经过训练的模型将未见过的数据投影到学习的嵌入空间中。

### 基本变换用法

```python
# Train on training data
trans = umap.UMAP(n_neighbors=15, random_state=42).fit(X_train)

# Transform test data
test_embedding = trans.transform(X_test)
```

### 与机器学习管道集成

```python
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import umap

# Split data
X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=0.2)

# Preprocess
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train UMAP
reducer = umap.UMAP(n_components=10, random_state=42)
X_train_embedded = reducer.fit_transform(X_train_scaled)
X_test_embedded = reducer.transform(X_test_scaled)

# Train classifier on embeddings
clf = SVC()
clf.fit(X_train_embedded, y_train)
accuracy = clf.score(X_test_embedded, y_test)
print(f"Test accuracy: {accuracy:.3f}")
```

### 重要考虑因素

**数据一致性：** 变换方法假设高维空间中的整体分布在训练数据和测试数据之间是一致的。当此假设失败时，请考虑改用参数化 UMAP。

**性能：** 转换操作非常高效（通常 <1 秒），但由于 Numba JIT 编译，初始调用可能会更慢。

**Scikit-learn 兼容性：** UMAP 遵循标准 sklearn 约定并在管道中无缝工作：

```python
from sklearn.pipeline import Pipeline

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('umap', umap.UMAP(n_components=10)),
    ('classifier', SVC())
])

pipeline.fit(X_train, y_train)
predictions = pipeline.predict(X_test)
```

## 高级功能

### 参数化 UMAP

参数化 UMAP 用学习的神经网络映射函数代替直接嵌入优化。

**与标准 UMAP 的主要区别：**
- 使用 TensorFlow/Keras 训练编码器网络
- 实现新数据的高效转换
- 支持通过解码器网络重建（逆变换）
- 允许自定义架构（用于图像的 CNN、用于序列的 RNN）

**安装：**
```bash
uv pip install umap-learn[parametric_umap]
# Requires TensorFlow 2.x
```

**基本用法：**
```python
from umap.parametric_umap import ParametricUMAP

# Default architecture (3-layer 100-neuron fully-connected network)
embedder = ParametricUMAP()
embedding = embedder.fit_transform(data)

# Transform new data efficiently
new_embedding = embedder.transform(new_data)
```

**定制架构：**
```python
import tensorflow as tf

# Define custom encoder
encoder = tf.keras.Sequential([
    tf.keras.layers.InputLayer(input_shape=(input_dim,)),
    tf.keras.layers.Dense(128, activation='relu'),
    tf.keras.layers.Dense(64, activation='relu'),
    tf.keras.layers.Dense(2)  # Output dimension
])

embedder = ParametricUMAP(encoder=encoder, dims=(input_dim,))
embedding = embedder.fit_transform(data)
```

**何时使用参数化 UMAP：**
- 训练后需要对新数据进行高效转换
- 需要重建能力（逆变换）
- 想要将 UMAP 与自动编码器结合起来
- 受益于专门架构，处理复杂的数据类型（图像、序列）

**何时使用标准 UMAP：**
- 需要简单和快速的原型设计
- 数据集很小，计算效率并不重要
- 不需要对未来数据进行学习转换

### 逆变换

逆变换可以从低维嵌入重建高维数据。

**基本用法：**
```python
reducer = umap.UMAP()
embedding = reducer.fit_transform(data)

# Reconstruct high-dimensional data from embedding coordinates
reconstructed = reducer.inverse_transform(embedding)
```

**重要限制：**
- 计算成本较高的操作
- 在嵌入的凸包之外效果不佳
- 集群之间存在差距的区域的准确度会降低

**使用案例：**
- 了解嵌入数据的结构
- 可视化集群之间的平滑过渡
- 探索数据点之间的插值
- 在嵌入空间中生成合成样本

**示例：探索嵌入空间：**
```python
import numpy as np

# Create grid of points in embedding space
x = np.linspace(embedding[:, 0].min(), embedding[:, 0].max(), 10)
y = np.linspace(embedding[:, 1].min(), embedding[:, 1].max(), 10)
xx, yy = np.meshgrid(x, y)
grid_points = np.c_[xx.ravel(), yy.ravel()]

# Reconstruct samples from grid
reconstructed_samples = reducer.inverse_transform(grid_points)
```

### 对齐UMAP

用于分析时间或相关数据集（例如时间序列实验、批量数据）：

```python
from umap import AlignedUMAP

# List of related datasets
datasets = [day1_data, day2_data, day3_data]

# Create aligned embeddings
mapper = AlignedUMAP().fit(datasets)
aligned_embeddings = mapper.embeddings_  # List of embeddings
```

**何时使用：** 比较相关数据集的嵌入，同时保持一致的坐标系。

## 再现性

为了确保结果可重现，请始终设置 `random_state` 参数：

```python
reducer = umap.UMAP(random_state=42)
```

UMAP 使用随机优化，因此在没有固定随机状态的情况下，运行之间的结果会略有不同。

## 常见问题及解决方案

**问题：** 组件断开连接或集群碎片
- **解决方案：** 增加 `n_neighbors` 以强调更多全局结构

**问题：** 簇过于分散或分离不佳
- **解决方案：** 减少 `min_dist` 以允许更紧密的包装

**问题：** 聚类结果不佳
- **解决方案：** 使用特定于聚类的参数（n_neighbors=30、min_dist=0.0、n_components=5-10）

**问题：** 转换结果与训练结果显着不同
- **解决方案：** 确保测试数据分布与训练相匹配，或使用参数化 UMAP

**问题：** 大型数据集上的性能缓慢
- **解决方案：** 设置`low_memory=True`（默认），或者先考虑用PCA降维

**问题：** 所有点都折叠成单个簇
- **解决方案：** 检查数据预处理（确保适当缩放），增加 `min_dist`

## 资源

###参考资料/

包含详细的API文档：
- `api_reference.md`：完整的UMAP类参数和方法

当需要详细参数信息或高级方法使用时，加载这些参考资料。