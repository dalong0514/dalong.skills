<!-- 此文件由机器翻译自 SKILL.md -->

---
姓名：永旺
描述：该技能应用于时间序列机器学习任务，包括分类、回归、聚类、预测、异常检测、分割和相似性搜索。当处理时态数据、顺序模式或时间索引观察需要标准 ML 方法之外的专用算法时使用。特别适合使用 scikit-learn 兼容 API 进行单变量和多变量时间序列分析。
---

# Aeon 时间序列机器学习

## 概述

Aeon 是一个与 scikit-learn 兼容的 Python 工具包，用于时间序列机器学习。它提供最先进的分类、回归、聚类、预测、异常检测、分割和相似性搜索算法。

## 何时使用此技能

在以下情况下应用此技能：
- 根据时间序列数据进行分类或预测
- 检测时间序列中的异常或变化点
- 聚类相似的时间序列模式
- 预测未来价值
- 寻找重复的模式（主题）或不寻常的子序列（不一致）
- 将时间序列与专门的距离度量进行比较
- 从时态数据中提取特征

## 安装

```bash
uv pip install aeon
```

## 核心能力

### 1.时间序列分类

将时间序列分类为预定义的类别。请参阅 `references/classification.md` 以获取完整的算法目录。

**快速入门：**
<<<代码块_1>>>

**算法选择：**
- **速度 + 性能**：`MiniRocketClassifier`、`Arsenal`
- **最大准确度**：`HIVECOTEV2`、`InceptionTimeClassifier`
- **可解释性**：`ShapeletTransformClassifier`、`Catch22Classifier`
- **小数据集**：`KNeighborsTimeSeriesClassifier` 具有 DTW 距离

### 2.时间序列回归

从时间序列预测连续值。请参阅 `references/regression.md` 了解算法。

**快速入门：**
<<<代码块_2>>>

### 3.时间序列聚类

对没有标签的相似时间序列进行分组。有关方法，请参阅 `references/clustering.md`。

**快速入门：**
<<<代码块_3>>>

### 4. 预测

预测未来的时间序列值。请参阅`references/forecasting.md`了解预报员。

**快速入门：**
<<<代码块_4>>>

### 5.异常检测

识别异常模式或异常值。有关检测器，请参阅 `references/anomaly_detection.md`。

**快速入门：**
<<<代码块_5>>>

### 6. 细分

将时间序列划分为具有变化点的区域。请参阅`references/segmentation.md`。

**快速入门：**
<<<代码块_6>>>

### 7. 相似性搜索

查找时间序列内或跨时间序列的相似模式。请参阅`references/similarity_search.md`。

**快速入门：**
```python
from aeon.similarity_search import StompMotif

# Find recurring patterns
motif_finder = StompMotif(window_size=50, k=3)
motifs = motif_finder.fit_predict(y)
```

## 特征提取和转换

转换时间序列以进行特征工程。请参阅`references/transformations.md`。

**火箭特点：**
```python
from aeon.transformations.collection.convolution_based import RocketTransformer

rocket = RocketTransformer()
X_features = rocket.fit_transform(X_train)

# Use features with any sklearn classifier
from sklearn.ensemble import RandomForestClassifier
clf = RandomForestClassifier()
clf.fit(X_features, y_train)
```

**统计特征：**
```python
from aeon.transformations.collection.feature_based import Catch22

catch22 = Catch22()
X_features = catch22.fit_transform(X_train)
```

**预处理：**
```python
from aeon.transformations.collection import MinMaxScaler, Normalizer

scaler = Normalizer()  # Z-normalization
X_normalized = scaler.fit_transform(X_train)
```

## 距离指标

专门的时间距离测量。请参阅 `references/distances.md` 了解完整目录。

**用途：**
```python
from aeon.distances import dtw_distance, dtw_pairwise_distance

# Single distance
distance = dtw_distance(x, y, window=0.1)

# Pairwise distances
distance_matrix = dtw_pairwise_distance(X_train)

# Use with classifiers
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier

clf = KNeighborsTimeSeriesClassifier(
    n_neighbors=5,
    distance="dtw",
    distance_params={"window": 0.2}
)
```

**可用距离：**
- **弹性**：DTW、DDTW、WDTW、ERP、EDR、LCSS、TWE、MSM
- **锁步**：欧几里得、曼哈顿、闵可夫斯基
- **基于形状**：形状 DTW、SBD

## 深度学习网络

时间序列的神经架构。请参阅`references/networks.md`。

**架构：**
- 卷积：`FCNClassifier`、`ResNetClassifier`、`InceptionTimeClassifier`
- 经常出现：`RecurrentNetwork`、`TCNNetwork`
- 自动编码器：`AEFCNClusterer`、`AEResNetClusterer`

**用途：**
```python
from aeon.classification.deep_learning import InceptionTimeClassifier

clf = InceptionTimeClassifier(n_epochs=100, batch_size=32)
clf.fit(X_train, y_train)
predictions = clf.predict(X_test)
```

## 数据集和基准测试

加载标准基准并评估性能。请参阅`references/datasets_benchmarking.md`。

**加载数据集：**
```python
from aeon.datasets import load_classification, load_regression

# Classification
X_train, y_train = load_classification("ArrowHead", split="train")

# Regression
X_train, y_train = load_regression("Covid3Month", split="train")
```

**基准测试：**
```python
from aeon.benchmarking import get_estimator_results

# Compare with published results
published = get_estimator_results("ROCKET", "GunPoint")
```

## 常见工作流程

### 分类管道

```python
from aeon.transformations.collection import Normalizer
from aeon.classification.convolution_based import RocketClassifier
from sklearn.pipeline import Pipeline

pipeline = Pipeline([
    ('normalize', Normalizer()),
    ('classify', RocketClassifier())
])

pipeline.fit(X_train, y_train)
accuracy = pipeline.score(X_test, y_test)
```

### 特征提取+传统机器学习

```python
from aeon.transformations.collection import RocketTransformer
from sklearn.ensemble import GradientBoostingClassifier

# Extract features
rocket = RocketTransformer()
X_train_features = rocket.fit_transform(X_train)
X_test_features = rocket.transform(X_test)

# Train traditional ML
clf = GradientBoostingClassifier()
clf.fit(X_train_features, y_train)
predictions = clf.predict(X_test_features)
```

### 可视化异常检测

```python
from aeon.anomaly_detection import STOMP
import matplotlib.pyplot as plt

detector = STOMP(window_size=50)
scores = detector.fit_predict(y)

plt.figure(figsize=(15, 5))
plt.subplot(2, 1, 1)
plt.plot(y, label='Time Series')
plt.subplot(2, 1, 2)
plt.plot(scores, label='Anomaly Scores', color='red')
plt.axhline(np.percentile(scores, 95), color='k', linestyle='--')
plt.show()
```

## 最佳实践

### 数据准备

1. **标准化**：大多数算法受益于 z 标准化
   ```python
   from aeon.transformations.collection import Normalizer
   normalizer = Normalizer()
   X_train = normalizer.fit_transform(X_train)
   X_test = normalizer.transform(X_test)
   ```

2. **处理缺失值**：分析前估算
   ```python
   from aeon.transformations.collection import SimpleImputer
   imputer = SimpleImputer(strategy='mean')
   X_train = imputer.fit_transform(X_train)
   ```

3. **检查数据格式**：Aeon 期望形状 `(n_samples, n_channels, n_timepoints)`

### 型号选择

1. **从简单开始**：在深度学习之前从 ROCKET 变体开始
2. **使用验证**：分割训练数据以进行超参数调优
3. **比较基线**：针对简单方法进行测试（1-NN 欧几里德、Naive）
4. **考虑资源**：ROCKET 提高速度，深度学习（如果 GPU 可用）

### 算法选择指南

**对于快速原型制作：**
- 分类：`MiniRocketClassifier`
- 回归：`MiniRocketRegressor`
- 聚类：`TimeSeriesKMeans` 与欧几里得

**为了获得最大准确度：**
- 分类：`HIVECOTEV2`、`InceptionTimeClassifier`
- 回归：`InceptionTimeRegressor`
- 预测：`ARIMA`、`TCNForecaster`

**为了可解释性：**
- 分类：`ShapeletTransformClassifier`、`Catch22Classifier`
- 功能：`Catch22`、`TSFresh`

**对于小数据集：**
- 基于距离：`KNeighborsTimeSeriesClassifier` 与 DTW
- 避免：深度学习（需要大数据）

## 参考文档

详细信息可参见`references/`：
- `classification.md` - 所有分类算法
- `regression.md` - 回归方法
- `clustering.md` - 聚类算法
- `forecasting.md` - 预测方法
- `anomaly_detection.md` - 异常检测方法
- `segmentation.md` - 分割算法
- `similarity_search.md` - 模式匹配和主题发现
- `transformations.md` - 特征提取和预处理
- `distances.md` - 时间序列距离度量
- `networks.md` - 深度学习架构
- `datasets_benchmarking.md` - 数据加载和评估工具

## 其他资源

- 文档：https://www.aeon-toolkit.org/
- GitHub：https://github.com/aeon-toolkit/aeon
- 示例：https://www.aeon-toolkit.org/en/stable/examples.html
- API 参考：https://www.aeon-toolkit.org/en/stable/api_reference.html