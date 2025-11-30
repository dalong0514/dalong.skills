<!-- 此文件由机器翻译自 anomaly_detection.md -->

# 异常检测

Aeon 提供异常检测方法，用于识别序列和集合级别的时间序列中的异常模式。

## 收集异常检测器

检测集合中的异常时间序列：

- `ClassificationAdapter` - 调整分类器以进行异常检测
  - 对正常数据进行训练，在预测期间标记异常值
  - **使用时**：已标记正常数据，需要基于分类的方法

- `OutlierDetectionAdapter` - 包裹 sklearn 离群值检测器
  - 适用于 IsolationForest、LOF、OneClassSVM
  - **使用时间**：想要在集合上使用 sklearn 异常检测器

## 系列异常探测器

检测单个时间序列内的异常点或子序列。

### 基于距离的方法

使用相似性度量来识别异常：

- `CBLOF` - 基于集群的局部异常值因子
  - 对数据进行聚类，根据聚类属性识别异常值
  - **使用时**：异常形成稀疏簇

- `KMeansAD` - 基于 K 均值的异常检测
  - 到最近聚类中心的距离表示异常
  - **使用时间**：正常模式聚类良好

- `LeftSTAMPi` - 左 STAMP 增量
  - 用于在线异常检测的矩阵配置文件
  - **使用时**：流数据，需要在线检测

- `STOMP` - 可扩展时间序列有序搜索矩阵配置文件
  - 计算子序列异常的矩阵轮廓
  - **使用时机**：Discord 发现、主题检测

- `MERLIN` - 基于矩阵配置文件的方法
  - 高效的矩阵轮廓计算
  - **使用时间**：大时间序列，需要可扩展性

- `LOF` - 适用于时间序列的局部异常值因子
  - 基于密度的异常值检测
  - **使用时**：低密度区域的异常

- `ROCKAD` - 基于 ROCKET 的半监督检测
  - 使用 ROCKET 特征进行异常识别
  - **使用时**：有一些标记数据，需要基于特征的方法

### 基于分布的方法

分析统计分布：

- `COPOD` - 基于 Copula 的异常值检测
  - 建立边际分布和联合分布模型
  - **使用时**：多维时间序列，复杂依赖关系

- `DWT_MLEAD` - 离散小波变换多级异常检测
  - 将系列分解为频段
  - **使用时间**：特定频率出现异常

### 基于隔离的方法

使用隔离原则：

- `IsolationForest` - 基于随机森林的隔离
  - 异常点比正常点更容易隔离
  - **使用时**：高维数据，没有关于分布的假设

- `OneClassSVM` - 用于新颖性检测的支持向量机
  - 学习正常数据周围的边界
  - **使用时**：明确定义的正常区域，需要稳健的边界

- `STRAY` - 流式稳健异常检测
  - 对数据分布变化具有鲁棒性
  - **使用时间**：流数据、分布变化

### 外部库集成

- `PyODAdapter` - 将 PyOD 库桥接到 aeon
  - 访问 40 多个 PyOD 异常检测器
  - **使用时**：需要特定的 PyOD 算法

## 快速入门

```python
from aeon.anomaly_detection import STOMP
import numpy as np

# Create time series with anomaly
y = np.concatenate([
    np.sin(np.linspace(0, 10, 100)),
    [5.0],  # Anomaly spike
    np.sin(np.linspace(10, 20, 100))
])

# Detect anomalies
detector = STOMP(window_size=10)
anomaly_scores = detector.fit_predict(y)

# Higher scores indicate more anomalous points
threshold = np.percentile(anomaly_scores, 95)
anomalies = anomaly_scores > threshold
```

## 点与子序列异常

- **点异常**：单个异常值
  - 使用：COPOD、DWT_MLEAD、IsolationForest

- **后续异常**（不一致）：不寻常的模式
  - 使用：STOMP、LeftSTAMPi、MERLIN

- **集体异常**：形成不寻常模式的点组
  - 使用：矩阵轮廓方法，基于聚类

## 评估指标

用于异常检测的专门指标：

<<<代码块_1>>>

## 算法选择

- **速度优先**：KMeansAD、IsolationForest
- **准确性优先**：STOMP、COPOD
- **流数据**：LeftSTAMPi、STRAY
- **不和谐发现**：STOMP、MERLIN
- **多维**：COPOD、PyODAdapter
- **半监督**：ROCKAD、OneClassSVM
- **无训练数据**：IsolationForest、STOMP

## 最佳实践

1. **标准化数据**：许多方法对尺度敏感
2. **选择窗口大小**：对于矩阵轮廓方法，窗口大小至关重要
3. **设置阈值**：使用基于百分位数或特定于域的阈值
4. **验证结果**：可视化检测以验证意义
5. **处理季节性**：在检测之前消除趋势/淡季化