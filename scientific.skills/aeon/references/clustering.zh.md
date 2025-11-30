<!-- 此文件由机器翻译自 clustering.md -->

# 时间序列聚类

Aeon 提供适用于时态数据的聚类算法，具有专门的距离度量和平均方法。

## 分区算法

适用于时间序列的标准 k-means/k-medoids：

- `TimeSeriesKMeans` - 具有时间距离度量的 K 均值（DTW、欧几里德等）
- `TimeSeriesKMedoids` - 使用实际时间序列作为聚类中心
- `TimeSeriesKShape` - 基于形状的聚类算法
- `TimeSeriesKernelKMeans` - 基于内核的非线性模式变体

**使用时**：已知簇数、预期的球形簇形状。

## 大数据集方法

大型集合的高效聚类：

- `TimeSeriesCLARA` - 通过采样对大型应用程序进行集群
- `TimeSeriesCLARANS` - CLARA 的随机搜索变体

**使用时间**：数据集对于标准 k-medoid 来说太大，需要可扩展性。

## 弹性距离聚类

专门用于基于对齐的相似性：

- `KASBA` - 具有平移不变弹性平均的 K 均值
- `ElasticSOM` - 使用弹性距离的自组织映射

**使用时间**：时间序列有时间变化或扭曲。

## 光谱方法

基于图的聚类：

- `KSpectralCentroid` - 具有质心计算的谱聚类

**使用时**：非凸簇形状，需要基于图的方法。

## 深度学习聚类

具有自动编码器的基于神经网络的聚类：

- `AEFCNClusterer` - 全卷积自动编码器
- `AEResNetClusterer` - 残差网络自动编码器
- `AEDCNNClusterer` - 扩张的 CNN 自动编码器
- `AEDRNNClusterer` - 扩张 RNN 自动编码器
- `AEBiGRUClusterer` - 双向 GRU 自动编码器
- `AEAttentionBiGRUClusterer` - 注意力增强型 BiGRU 自动编码器

**使用时机**：大型数据集、需要学习的表示或复杂的模式。

## 基于特征的聚类

聚类前变换到特征空间：

- `Catch22Clusterer` - 22 个规范特征上的聚类
- `SummaryClusterer` - 使用汇总统计信息
- `TSFreshClusterer` - 自动 tsfresh 功能

**使用时间**：原始时间序列信息不丰富，需要可解释的功能。

## 成分

构建自定义集群管道：

- `ClustererPipeline` - 带集群器的链式变压器

## 平均方法

计算时间序列的聚类中心：

- `mean_average` - 算术平均值
- `ba_average` - 使用 DTW 进行重心平均
- `kasba_average` - 平移不变平均
- `shift_invariant_average` - 一般移位不变方法

**使用时**：需要代表性聚类中心进行可视化或初始化。

## 快速入门

```python
from aeon.clustering import TimeSeriesKMeans
from aeon.datasets import load_classification

# Load data (using classification data for clustering)
X_train, _ = load_classification("GunPoint", split="train")

# Cluster time series
clusterer = TimeSeriesKMeans(
    n_clusters=3,
    distance="dtw",  # Use DTW distance
    averaging_method="ba"  # Barycentric averaging
)
labels = clusterer.fit_predict(X_train)
centers = clusterer.cluster_centers_
```

## 算法选择

- **速度优先**：TimeSeriesKMeans 与欧几里德距离
- **时间对齐**：KASBA、TimeSeriesKMeans 和 DTW
- **大型数据集**：TimeSeriesCLARA、TimeSeriesCLARANS
- **复杂模式**：深度学习聚类器
- **可解释性**：Catch22Clusterer、SummaryClusterer
- **非凸簇**：KSpectralCentroid

## 距离指标

兼容的距离指标包括：
- 欧几里得、曼哈顿、闵可夫斯基（锁步）
- DTW、DDTW、WDTW（弹性对齐）
- ERP、EDR、LCSS（基于编辑）
- MSM、TWE（特种弹力）

## 评价

使用 sklearn 或 aeon 基准测试中的聚类指标：
- 剪影得分
- 戴维斯-布尔丁指数
- 卡林斯基-哈拉巴斯指数