<!-- 此文件由机器翻译自 unsupervised_learning.md -->

# 无监督学习参考

## 概述

无监督学习通过聚类、降维和密度估计来发现未标记数据中的模式。

## 聚类

### K 均值

**KMeans (`sklearn.cluster.KMeans`)**
- 基于分区的聚类为 K 个簇
- 关键参数：
  - `n_clusters`：要形成的簇数
  - `init`：初始化方法（'k-means++'、'随机'）
  - `n_init`：初始化次数（默认=10）
  - `max_iter`：最大迭代次数
- 使用场合：了解簇数、球形簇形状
- 快速且可扩展
- 示例：
```python
from sklearn.cluster import KMeans

model = KMeans(n_clusters=3, init='k-means++', n_init=10, random_state=42)
labels = model.fit_predict(X)
centers = model.cluster_centers_

# Inertia (sum of squared distances to nearest center)
print(f"Inertia: {model.inertia_}")
```

**MiniBatchKMeans**
- 使用小批量更快的 K-Means
- 使用时机：大型数据集，需要更快的训练
- 比 K-Means 准确度稍低
- 示例：
<<<代码块_1>>>

### 基于密度的聚类

**DBSCAN (`sklearn.cluster.DBSCAN`)**
- 基于密度的空间聚类
- 关键参数：
  - `eps`：两个相邻样本之间的最大距离
  - `min_samples`：邻域内形成核心点的最小样本
  - `metric`：距离度量
- 在以下情况下使用：任意簇形状、存在噪声/异常值
- 自动确定簇的数量
- 将噪声点标记为-1
- 示例：
<<<代码块_2>>>

**HDBSCAN (`sklearn.cluster.HDBSCAN`)**
- 具有自适应 epsilon 的分层 DBSCAN
- 比 DBSCAN 更稳健
- 关键参数：`min_cluster_size`
- 使用场合：不同密度的簇
- 示例：
<<<代码块_3>>>

**光学 (`sklearn.cluster.OPTICS`)**
- 排序点以识别聚类结构
- 与 DBSCAN 类似，但不需要 eps 参数
- 关键参数：`min_samples`、`max_eps`
- 使用场合：变化密度、探索性分析
- 示例：
<<<代码块_4>>>

### 层次聚类

**凝聚聚类**
- 自下而上的层次聚类
- 关键参数：
  - `n_clusters`：簇数（或使用`distance_threshold`）
  - `linkage`: '病房', '完整', '平均', '单一'
  - `metric`：距离度量
- 使用时：需要树状图，层次结构很重要
- 示例：
<<<代码块_5>>>

### 其他聚类方法

**均值平移**
- 通过将点移向密度模式来查找簇
- 自动确定簇的数量
- 关键参数：`bandwidth`
- 使用时：不知道簇的数量、任意形状
- 示例：
<<<代码块_6>>>

**光谱聚类**
- 使用基于图的方法和特征值
- 关键参数：`n_clusters`、`affinity`（'rbf'、'nearest_neighbors'）
- 使用场合：非凸簇、图结构
- 示例：
```python
from sklearn.cluster import SpectralClustering

model = SpectralClustering(n_clusters=3, affinity='rbf', random_state=42)
labels = model.fit_predict(X)
```

**亲和力传播**
- 通过消息传递查找范例
- 自动确定簇的数量
- 关键参数：`damping`、`preference`
- 使用时：不知道簇的数量
- 示例：
```python
from sklearn.cluster import AffinityPropagation

model = AffinityPropagation(damping=0.9, random_state=42)
labels = model.fit_predict(X)
n_clusters = len(model.cluster_centers_indices_)
```

**桦木**
- 使用层次结构平衡迭代减少和聚类
- 大型数据集的内存效率
- 关键参数：`n_clusters`、`threshold`、`branching_factor`
- 使用时机：非常大的数据集
- 示例：
```python
from sklearn.cluster import Birch

model = Birch(n_clusters=3, threshold=0.5)
labels = model.fit_predict(X)
```

### 聚类评估

**已知基本事实时的指标：**
```python
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics import adjusted_mutual_info_score, fowlkes_mallows_score

# Compare predicted labels with true labels
ari = adjusted_rand_score(y_true, y_pred)
nmi = normalized_mutual_info_score(y_true, y_pred)
ami = adjusted_mutual_info_score(y_true, y_pred)
fmi = fowlkes_mallows_score(y_true, y_pred)
```

**没有基本事实的指标：**
```python
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.metrics import davies_bouldin_score

# Silhouette: [-1, 1], higher is better
silhouette = silhouette_score(X, labels)

# Calinski-Harabasz: higher is better
ch_score = calinski_harabasz_score(X, labels)

# Davies-Bouldin: lower is better
db_score = davies_bouldin_score(X, labels)
```

**K-Means 的肘部法：**
```python
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

inertias = []
K_range = range(2, 11)
for k in K_range:
    model = KMeans(n_clusters=k, random_state=42)
    model.fit(X)
    inertias.append(model.inertia_)

plt.plot(K_range, inertias, 'bo-')
plt.xlabel('Number of clusters')
plt.ylabel('Inertia')
plt.title('Elbow Method')
```

## 降维

### 主成分分析 (PCA)

**PCA (`sklearn.decomposition.PCA`)**
- 使用 SVD 进行线性降维
- 关键参数：
  - `n_components`：分量数量（解释方差的 int 或 float）
  - `whiten`：将组件白化为单位方差
- 使用时：线性关系，想要解释方差
- 示例：
```python
from sklearn.decomposition import PCA

# Keep components explaining 95% variance
pca = PCA(n_components=0.95)
X_reduced = pca.fit_transform(X)

print(f"Original dimensions: {X.shape[1]}")
print(f"Reduced dimensions: {X_reduced.shape[1]}")
print(f"Explained variance ratio: {pca.explained_variance_ratio_}")
print(f"Total variance explained: {pca.explained_variance_ratio_.sum()}")

# Or specify exact number of components
pca = PCA(n_components=2)
X_2d = pca.fit_transform(X)
```

**增量PCA**
- PCA 适用于无法放入内存的大型数据集
- 批量处理数据
- 关键参数：`n_components`、`batch_size`
- 示例：
```python
from sklearn.decomposition import IncrementalPCA

pca = IncrementalPCA(n_components=50, batch_size=100)
X_reduced = pca.fit_transform(X)
```

**内核PCA**
- 使用内核进行非线性降维
- 关键参数：`n_components`、`kernel`（'线性'、'poly'、'rbf'、'sigmoid'）
- 使用时：非线性关系
- 示例：
```python
from sklearn.decomposition import KernelPCA

pca = KernelPCA(n_components=2, kernel='rbf', gamma=0.1)
X_reduced = pca.fit_transform(X)
```

### 流形学习

**t-SNE (`sklearn.manifold.TSNE`)**
- t-分布随机邻域嵌入
- 非常适合 2D/3D 可视化
- 关键参数：
  - `n_components`：通常为 2 或 3
  - `perplexity`：局部和全局结构之间的平衡（5-50）
- `learning_rate`：通常为 10-1000
  - `n_iter`：迭代次数（最少 250）
- 使用场合：可视化高维数据
- 注意：在大型数据集上速度很慢，没有transform()方法
- 示例：
```python
from sklearn.manifold import TSNE

tsne = TSNE(n_components=2, perplexity=30, learning_rate=200, n_iter=1000, random_state=42)
X_embedded = tsne.fit_transform(X)

# Visualize
import matplotlib.pyplot as plt
plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=labels, cmap='viridis')
plt.title('t-SNE visualization')
```

**UMAP（不在 scikit-learn 中，但兼容）**
- 均匀流形逼近和投影
- 比 t-SNE 更快，更好地保留全局结构
- 安装：`uv pip install umap-learn`
- 示例：
```python
from umap import UMAP

reducer = UMAP(n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
X_embedded = reducer.fit_transform(X)
```

**等距图**
- 等轴测图
- 保留测地距离
- 关键参数：`n_components`、`n_neighbors`
- 使用场合：非线性流形
- 示例：
```python
from sklearn.manifold import Isomap

isomap = Isomap(n_components=2, n_neighbors=5)
X_embedded = isomap.fit_transform(X)
```

**局部线性嵌入 (LLE)**
- 保留当地邻里结构
- 关键参数：`n_components`、`n_neighbors`
- 示例：
```python
from sklearn.manifold import LocallyLinearEmbedding

lle = LocallyLinearEmbedding(n_components=2, n_neighbors=10)
X_embedded = lle.fit_transform(X)
```

**MDS（多维尺度）**
- 保留成对距离
- 关键参数：`n_components`、`metric`（真/假）
- 示例：
```python
from sklearn.manifold import MDS

mds = MDS(n_components=2, metric=True, random_state=42)
X_embedded = mds.fit_transform(X)
```

### 矩阵分解

**NMF（非负矩阵分解）**
- 分解为非负矩阵
- 关键参数：`n_components`、`init`（'nndsvd'、'随机'）
- 使用时：数据非负（图像、文本）
- 可解释的组件
- 示例：
```python
from sklearn.decomposition import NMF

nmf = NMF(n_components=10, init='nndsvd', random_state=42)
W = nmf.fit_transform(X)  # Document-topic matrix
H = nmf.components_  # Topic-word matrix
```

**截断SVD**
- 稀疏矩阵的 SVD
- 与 PCA 类似，但适用于稀疏数据
- 使用场合：文本数据、稀疏矩阵
- 示例：
```python
from sklearn.decomposition import TruncatedSVD

svd = TruncatedSVD(n_components=100, random_state=42)
X_reduced = svd.fit_transform(X_sparse)
print(f"Explained variance: {svd.explained_variance_ratio_.sum()}")
```

**快速ICA**
- 独立成分分析
- 将多变量信号分离成独立的分量
- 关键参数：`n_components`
- 使用场合：信号分离（例如音频、EEG）
- 示例：
```python
from sklearn.decomposition import FastICA

ica = FastICA(n_components=10, random_state=42)
S = ica.fit_transform(X)  # Independent sources
A = ica.mixing_  # Mixing matrix
```

**潜在狄利克雷分配 (LDA)**
- 文本数据的主题建模
- 关键参数：`n_components`（主题数量）、`learning_method`（'批量'、'在线'）
- 使用场合：主题建模、文档聚类
- 示例：
```python
from sklearn.decomposition import LatentDirichletAllocation

lda = LatentDirichletAllocation(n_components=10, random_state=42)
doc_topics = lda.fit_transform(X_counts)  # Document-topic distribution

# Get top words for each topic
feature_names = vectorizer.get_feature_names_out()
for topic_idx, topic in enumerate(lda.components_):
    top_words = [feature_names[i] for i in topic.argsort()[-10:]]
    print(f"Topic {topic_idx}: {', '.join(top_words)}")
```

## 异常值和新颖性检测

### 异常值检测

**隔离森林**
- 使用随机树隔离异常
- 关键参数：
  - `contamination`：异常值的预期比例
  - `n_estimators`：树的数量
- 使用场合：高维数据，效率重要
- 示例：
```python
from sklearn.ensemble import IsolationForest

model = IsolationForest(contamination=0.1, random_state=42)
predictions = model.fit_predict(X)  # -1 for outliers, 1 for inliers
```

**局部异常值因素**
- 测量局部密度偏差
- 关键参数：`n_neighbors`、`contamination`
- 使用场合：不同密度区域
- 示例：
```python
from sklearn.neighbors import LocalOutlierFactor

lof = LocalOutlierFactor(n_neighbors=20, contamination=0.1)
predictions = lof.fit_predict(X)  # -1 for outliers, 1 for inliers
outlier_scores = lof.negative_outlier_factor_
```

**一类 SVM**
- 学习正常数据周围的决策边界
- 关键参数：`nu`（离群值上限）、`kernel`、`gamma`
- 使用场合：正常数据的小型训练集
- 示例：
```python
from sklearn.svm import OneClassSVM

model = OneClassSVM(nu=0.1, kernel='rbf', gamma='auto')
model.fit(X_train)
predictions = model.predict(X_test)  # -1 for outliers, 1 for inliers
```

**椭圆包络**
- 假设高斯分布
- 关键参数：`contamination`
- 使用时：数据呈高斯分布
- 示例：
```python
from sklearn.covariance import EllipticEnvelope

model = EllipticEnvelope(contamination=0.1, random_state=42)
predictions = model.fit_predict(X)
```

## 高斯混合模型

**高斯混合**
- 混合高斯概率聚类
- 关键参数：
  - `n_components`：混合成分的数量
  - `covariance_type`：“完整”、“捆绑”、“诊断”、“球形”
- 使用时：软聚类，需要概率估计
- 示例：
```python
from sklearn.mixture import GaussianMixture

gmm = GaussianMixture(n_components=3, covariance_type='full', random_state=42)
gmm.fit(X)

# Predict cluster labels
labels = gmm.predict(X)

# Get probability of each cluster
probabilities = gmm.predict_proba(X)

# Information criteria for model selection
print(f"BIC: {gmm.bic(X)}")  # Lower is better
print(f"AIC: {gmm.aic(X)}")  # Lower is better
```

## 选择正确的方法

### 聚类：
- **了解 K，球形簇**：K 均值
- **任意形状、噪声**：DBSCAN、HDBSCAN
- **层次结构**：凝聚聚类
- **非常大的数据**：MiniBatchKMeans、BIRCH
- **概率**：GaussianMixture

### 降维：
- **线性、方差解释**：PCA
- **非线性、可视化**：t-SNE、UMAP
- **非负数据**：NMF
- **稀疏数据**：TruncatedSVD
- **主题建模**：LatentDirichletAllocation

### 异常值检测：
- **高维**：IsolationForest
- **不同密度**：LocalOutlierFactor
- **高斯数据**：EllipticEnvelope