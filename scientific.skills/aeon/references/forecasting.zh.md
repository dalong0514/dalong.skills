<!-- 此文件由机器翻译自 forecasting.md -->

# 时间序列预测

Aeon 提供用于预测未来时间序列值的预测算法。

## 简单方法和基线方法

用于比较的简单预测策略：

- `NaiveForecaster` - 多种策略：最后值、平均值、季节天真
  - 参数：`strategy`（“最后”、“平均值”、“季节性”）、`sp`（季节性时段）
  - **使用时间**：建立基线或简单模式

## 统计模型

经典时间序列预测方法：

### ARIMA
- `ARIMA` - 自回归综合移动平均线
  - 参数：`p`（AR 顺序）、`d`（差分）、`q`（MA 顺序）
  - **使用时**：线性模式、平稳或差分平稳序列

### 指数平滑
- `ETS` - 错误-趋势-季节分解
  - 参数：`error`、`trend`、`seasonal` 类型
  - **使用时间**：存在趋势和季节性模式

### 阈值自回归
- `TAR` - 用于状态切换的阈值自回归模型
- `AutoTAR` - 自动阈值发现
  - **使用时间**：系列在不同的情况下表现出不同的行为

### 西塔法
- `Theta` - 经典 Theta 预测
  - 参数：`theta`、`weights`用于分解
  - **使用时间**：需要简单但有效的基线

### 时变参数
- `TVP` - 带卡尔曼滤波的时变参数模型
  - **使用时间**：参数随时间变化

## 深度学习预测器

用于复杂时间模式的神经网络：

- `TCNForecaster` - 时间卷积网络
  - 大感受野的扩张卷积
  - **使用时间**：长序列，需要非循环架构

- `DeepARNetwork` - 使用 RNN 进行概率预测
  - 提供预测区间
  - **使用时**：需要概率预测、不确定性量化

## 基于回归的预测

将回归应用于滞后特征：

- `RegressionForecaster` - 包装回归量以进行预测
  - 参数：`window_length`、`horizon`
  - **使用时间**：想要使用任何回归器作为预测器

## 快速入门

```python
from aeon.forecasting.naive import NaiveForecaster
from aeon.forecasting.arima import ARIMA
import numpy as np

# Create time series
y = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

# Naive baseline
naive = NaiveForecaster(strategy="last")
naive.fit(y)
forecast_naive = naive.predict(fh=[1, 2, 3])

# ARIMA model
arima = ARIMA(order=(1, 1, 1))
arima.fit(y)
forecast_arima = arima.predict(fh=[1, 2, 3])
```

## 预测范围

预测范围 (`fh`) 指定要预测的未来时间点：

<<<代码块_1>>>

## 型号选择

- **基线**：具有季节性策略的 NaiveForecaster
- **线性模式**：ARIMA
- **趋势 + 季节性**：ETS
- **制度变化**：TAR、AutoTAR
- **复杂模式**：TCNForecaster
- **概率**：DeepARNetwork
- **长序列**：TCNForecaster
- **短序列**：ARIMA、ETS

## 评估指标

使用标准预测指标：

<<<代码块_2>>>

## 外生变量

许多预测者支持外生特征：

<<<代码块_3>>>

## 基类

- `BaseForecaster` - 所有预测者的抽象基础
- `BaseDeepForecaster` - 深度学习预测器的基础

扩展它们以实现自定义预测算法。