<!-- 此文件由机器翻译自 segmentation.md -->

# 时间序列分割

Aeon 提供算法将时间序列划分为具有不同特征的区域，识别变化点和边界。

## 分割算法

### 二值分割
- `BinSegmenter` - 递归二进制分段
  - 在最重要的变化点迭代地分割系列
  - 参数：`n_segments`、`cost_function`
  - **使用时**：已知段数，分层结构

### 基于分类
- `ClaSPSegmenter` - 分类分数配置文件
  - 使用分类性能来识别边界
  - 发现分类区分邻居的分段
  - **使用时间**：分段具有不同的时间模式

### 基于快速模式
- `FLUSSSegmenter` - 快速低成本单能语义分割
  - 使用弧交叉进行高效语义分割
  - 基于矩阵轮廓
  - **使用时间**：大时间序列，需要速度和模式发现

### 信息论
- `InformationGainSegmenter` - 信息增益最大化
  - 找到最大化信息增益的边界
  - **使用时间**：段之间的统计差异

### 高斯建模
- `GreedyGaussianSegmenter` - 贪婪高斯近似
  - 将分段建模为高斯分布
  - 增量添加变化点
  - **使用时**：段遵循高斯分布

### 分层聚合
- `EAggloSegmenter` - 自下而上的合并方法
  - 通过聚集估计变化点
  - **使用时**：想要分层分段结构

### 隐马尔可夫模型
- `HMMSegmenter` - 带维特比解码的 HMM
  - 基于概率状态的分割
  - **使用时**：段代表隐藏状态

### 基于维度
- `HidalgoSegmenter` - 异构内在维数算法
  - 检测局部维度的变化
  - **使用时间**：段之间的维度变化

### 基线
- `RandomSegmenter` - 随机变化点生成
  - **使用时**：需要原假设基线

## 快速入门

```python
from aeon.segmentation import ClaSPSegmenter
import numpy as np

# Create time series with regime changes
y = np.concatenate([
    np.sin(np.linspace(0, 10, 100)),      # Segment 1
    np.cos(np.linspace(0, 10, 100)),      # Segment 2
    np.sin(2 * np.linspace(0, 10, 100))   # Segment 3
])

# Segment the series
segmenter = ClaSPSegmenter()
change_points = segmenter.fit_predict(y)

print(f"Detected change points: {change_points}")
```

## 输出格式

分段器返回变化点索引：

<<<代码块_1>>>

## 算法选择

- **速度优先**：FLUSSSegmenter、BinSegmenter
- **精度优先**：ClaSPSegmenter、HMMSegmenter
- **已知段计数**：带有 n_segments 参数的 BinSegmenter
- **未知段计数**：ClaSPSegmenter、InformationGainSegmenter
- **模式更改**：FLUSSSegmenter、ClaSPSegmenter
- **统计变化**：InformationGainSegmenter、GreedyGaussianSegmenter
- **状态转换**：HMMSegmenter

## 常见用例

### 政权变更检测
确定时间序列行为何时发生根本性变化：

<<<代码块_2>>>

### 活动细分
将传感器数据细分为活动：

<<<代码块_3>>>

### 季节性边界检测
查找时间序列中的季节转换：

<<<代码块_4>>>

## 评估指标

使用分段质量指标：

<<<代码块_5>>>

## 最佳实践

1. **标准化数据**：确保变化检测不受规模支配
2. **选择合适的指标**：不同的算法优化不同的标准
3. **验证段**：可视化以验证有意义的边界
4. **处理噪声**：在分割之前考虑平滑
5. **领域知识**：使用预期的段计数（如果已知）
6. **参数调整**：调整灵敏度参数（阈值、惩罚）

## 可视化

<<<代码块_6>>>