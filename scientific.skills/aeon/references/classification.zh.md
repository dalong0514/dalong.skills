<!-- 此文件由机器翻译自 classification.md -->

# 时间序列分类

Aeon 提供 13 类时间序列分类器以及 scikit-learn 兼容的 API。

## 基于卷积的分类器

应用随机卷积变换进行有效的特征提取：

- `Arsenal` - 具有不同内核的 ROCKET 分类器的集合
- `HydraClassifier` - 带扩张的多分辨率卷积
- `RocketClassifier` - 具有岭回归的随机卷积核
- `MiniRocketClassifier` - 简化的 ROCKET 变体以提高速度
- `MultiRocketClassifier` - 组合多个 ROCKET 变体

**使用时机**：需要快速、可扩展的分类，并且在不同的数据集上具有强大的性能。

## 深度学习分类器

针对时间序列优化的神经网络架构：

- `FCNClassifier` - 全卷积网络
- `ResNetClassifier` - 具有跳过连接的剩余网络
- `InceptionTimeClassifier` - 多尺度初始模块
- `TimeCNNClassifier` - 时间序列的标准 CNN
- `MLPClassifier` - 多层感知器基线
- `EncoderClassifier` - 通用编码器包装器
- `DisjointCNNClassifier` - 以 Shapelet 为中心的架构

**使用时机**：可用的大型数据集、需要端到端学习或复杂的时间模式。

## 基于字典的分类器

将时间序列转换为符号表示：

- `BOSSEnsemble` - 带整体投票的 Bag-of-SFA-Symbols
- `TemporalDictionaryEnsemble` - 多种字典方法组合
- `WEASEL` - 用于时间序列分类的字提取操作
- `MrSEQLClassifier` - 多重符号序列学习

**何时使用**：需要可解释的模型、稀疏模式或符号推理。

## 基于距离的分类器

利用专门的时间序列距离度量：

- `KNeighborsTimeSeriesClassifier` - 具有时间距离的 k-NN（DTW、LCSS、ERP 等）
- `ElasticEnsemble` - 组合多个弹性距离测量
- `ProximityForest` - 使用基于距离的分割的树集成

**使用时机**：数据集较小，需要基于相似性的分类或可解释的决策。

## 基于特征的分类器

分类前提取统计和签名特征：

- `Catch22Classifier` - 22 个规范时间序列特征
- `TSFreshClassifier` - 通过 tsfresh 自动提取特征
- `SignatureClassifier` - 路径签名转换
- `SummaryClassifier` - 摘要统计提取
- `FreshPRINCEClassifier` - 组合多个特征提取器

**何时使用**：需要可解释的功能、可用的领域专业知识或功能工程方法。

## 基于区间的分类器

从随机或监督间隔中提取特征：

- `CanonicalIntervalForestClassifier` - 带有决策树的随机区间特征
- `DrCIFClassifier` - 具有 catch22 功能的多样化表示 CIF
- `TimeSeriesForestClassifier` - 具有汇总统计信息的随机间隔
- `RandomIntervalClassifier` - 简单的基于间隔的方法
- `RandomIntervalSpectralEnsembleClassifier` - 区间的频谱特征
- `SupervisedTimeSeriesForest` - 监督间隔选择

**使用时间**：区分模式发生在特定时间窗口内。

## 基于 Shapelet 的分类器

识别判别子序列（shapelet）：

- `ShapeletTransformClassifier` - 发现并使用有区别的 shapelet
- `LearningShapeletClassifier` - 通过梯度下降学习 shapelet
- `SASTClassifier` - 可扩展的近似 shapelet 变换
- `RDSTClassifier` - 随机扩张 shapelet 变换

**使用时**：需要可解释的判别模式或相位不变特征。

## 混合分类器

结合多种分类范式：

- `HIVECOTEV1` - 基于转换的集成的分层投票集体（版本 1）
- `HIVECOTEV2` - 具有更新组件的增强版本

**使用时间**：需要最大精度，可用计算资源。

## 早期分类

在观察整个时间序列之前进行预测：

- `TEASER` - 两层早期准确系列分类器
- `ProbabilityThresholdEarlyClassifier` - 置信度超过阈值时的预测

**使用时间**：需要实时决策，或者观察需要成本。

## 序数分类

处理有序的类标签：

- `OrdinalTDE` - 用于序数输出的时态字典集合

**使用时间**：类具有自然顺序（例如，严重性级别）。

## 合成工具

构建自定义管道和集成：

- `ClassifierPipeline` - 带分类器的链式变压器
- `WeightedEnsembleClassifier` - 分类器的加权组合
- `SklearnClassifierWrapper` - 调整 sklearn 分类器以适应时间序列

## 快速入门

```python
from aeon.classification.convolution_based import RocketClassifier
from aeon.datasets import load_classification

# Load data
X_train, y_train = load_classification("GunPoint", split="train")
X_test, y_test = load_classification("GunPoint", split="test")

# Train and predict
clf = RocketClassifier()
clf.fit(X_train, y_train)
accuracy = clf.score(X_test, y_test)
```

## 算法选择

- **速度优先**：MiniRocketClassifier，Arsenal
- **准确性优先**：HIVECOTEV2、InceptionTimeClassifier
- **可解释性**：ShapeletTransformClassifier、Catch22Classifier
- **小数据**：KNeighborsTimeSeriesClassifier，基于距离的方法
- **大数据**：深度学习分类器、ROCKET 变体