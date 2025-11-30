<!-- 此文件由机器翻译自 regression.md -->

# 时间序列回归

Aeon 提供跨 9 个类别的时间序列回归器，用于从时间序列预测连续值。

## 基于卷积的回归器

应用卷积核进行特征提取：

- `HydraRegressor` - 多分辨率扩张卷积
- `RocketRegressor` - 随机卷积核
- `MiniRocketRegressor` - 简化 ROCKET 以提高速度
- `MultiRocketRegressor` - 组合火箭变体
- `MultiRocketHydraRegressor` - 合并 ROCKET 和 Hydra 方法

**何时使用**：需要具有强大基线性能的快速回归。

## 深度学习回归器

用于端到端时间回归的神经架构：

- `FCNRegressor` - 全卷积网络
- `ResNetRegressor` - 具有跳过连接的剩余块
- `InceptionTimeRegressor` - 多尺度初始模块
- `TimeCNNRegressor` - 标准 CNN 架构
- `RecurrentRegressor` - RNN/LSTM/GRU 变体
- `MLPRegressor` - 多层感知器
- `EncoderRegressor` - 通用编码器包装器
- `LITERegressor` - 轻量级初始时间集合
- `DisjointCNNRegressor` - 专门的 CNN 架构

**使用时机**：大型数据集、复杂模式或需要特征学习。

## 基于距离的回归器

具有时间距离度量的 k 最近邻：

- `KNeighborsTimeSeriesRegressor` - 具有 DTW、LCSS、ERP 或其他距离的 k-NN

**使用时机**：小型数据集、局部相似性模式或可解释的预测。

## 基于特征的回归器

回归前提取统计特征：

- `Catch22Regressor` - 22 个规范时间序列特征
- `FreshPRINCERegressor` - 组合多个特征提取器的管道
- `SummaryRegressor` - 摘要统计功能
- `TSFreshRegressor` - 自动 tsfresh 特征提取

**何时使用**：需要可解释的功能或特定领域的功能工程。

## 混合回归器

结合多种方法：

- `RISTRegressor` - 随机区间 Shapelet 变换

**使用时间**：从组合间隔和 shapelet 方法中受益。

## 基于区间的回归器

从时间间隔中提取特征：

- `CanonicalIntervalForestRegressor` - 决策树的随机间隔
- `DrCIFRegressor` - 多样化表示 CIF
- `TimeSeriesForestRegressor` - 随机间隔集合
- `RandomIntervalRegressor` - 简单的基于间隔的方法
- `RandomIntervalSpectralEnsembleRegressor` - 光谱间隔特征
- `QUANTRegressor` - 基于分位数的区间特征

**使用时间**：预测模式发生在特定时间窗口内。

## 基于 Shapelet 的回归器

使用判别子序列进行预测：

- `RDSTRegressor` - 随机扩张 Shapelet 变换

**使用时**：需要相位不变的判别模式。

## 合成工具

构建自定义回归管道：

- `RegressorPipeline` - 带回归器的链式变压器
- `RegressorEnsemble` - 具有可学习权重的加权集成
- `SklearnRegressorWrapper` - 调整 sklearn 回归器以适应时间序列

## 实用程序

- `DummyRegressor` - 基线策略（平均值、中位数）
- `BaseRegressor` - 自定义回归器的抽象基础
- `BaseDeepRegressor` - 深度学习回归器的基础

## 快速入门

```python
from aeon.regression.convolution_based import RocketRegressor
from aeon.datasets import load_regression

# Load data
X_train, y_train = load_regression("Covid3Month", split="train")
X_test, y_test = load_regression("Covid3Month", split="test")

# Train and predict
reg = RocketRegressor()
reg.fit(X_train, y_train)
predictions = reg.predict(X_test)
```

## 算法选择

- **速度优先**：MiniRocketRegressor
- **精度优先**：InceptionTimeRegressor、MultiRocketHydraRegressor
- **可解释性**：Catch22Regressor、SummaryRegressor
- **小数据**：KNeighborsTimeSeriesRegressor
- **大数据**：深度学习回归器、ROCKET 变体
- **间隔模式**：DrCIFRegressor、CanonicalIntervalForestRegressor