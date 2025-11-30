<!-- 此文件由机器翻译自 similarity_search.md -->

# 相似度搜索

Aeon 提供了在时间序列内和跨时间序列查找相似模式的工具，包括子序列搜索、主题发现和近似最近邻。

## 子序列最近邻 (SNN)

查找时间序列中最相似的子序列。

### MASS 算法
- `MassSNN` - Mueen 相似性搜索算法
  - 快速归一化互相关以实现相似性
  - 有效计算距离剖面
  - **使用时**：需要精确的最近邻距离，大系列

### 基于 STOMP 的基序发现
- `StompMotif` - 发现重复出现的模式（主题）
  - 查找前k个最相似的子序列对
  - 基于矩阵轮廓计算
  - **使用时**：想要发现重复的模式

### 暴力破解基线
- `DummySNN` - 详尽的距离计算
  - 计算所有成对距离
  - **使用时**：小系列，需要精确的基线

## 集合级搜索

跨集合查找相似的时间序列。

### 近似最近邻（ANN）
- `RandomProjectionIndexANN` - 位置敏感哈希
  - 使用具有余弦相似度的随机投影
  - 建立索引以进行快速近似搜索
  - **使用时**：大量收集，速度比准确性更重要

## 快速入门：基序发现

```python
from aeon.similarity_search import StompMotif
import numpy as np

# Create time series with repeated patterns
pattern = np.sin(np.linspace(0, 2*np.pi, 50))
y = np.concatenate([
    pattern + np.random.normal(0, 0.1, 50),
    np.random.normal(0, 1, 100),
    pattern + np.random.normal(0, 0.1, 50),
    np.random.normal(0, 1, 100)
])

# Find top-3 motifs
motif_finder = StompMotif(window_size=50, k=3)
motifs = motif_finder.fit_predict(y)

# motifs contains indices of motif occurrences
for i, (idx1, idx2) in enumerate(motifs):
    print(f"Motif {i+1} at positions {idx1} and {idx2}")
```

## 快速入门：子序列搜索

<<<代码块_1>>>

## 快速入门：集合上的近似 NN

<<<代码块_2>>>

## 矩阵简介

矩阵轮廓是许多相似性搜索任务的基本数据结构：

- **距离概况**：从查询到所有子序列的距离
- **矩阵轮廓**：每个子序列到任何其他子序列的最小距离
- **Motif**：具有最小距离的子序列对
- **不和谐**：具有最大最小距离的子序列（异常）

<<<代码块_3>>>

## 算法选择

- **精确子序列搜索**：MassSNN
- **主题发现**：StompMotif
- **异常检测**：矩阵配置文件（参见 anomaly_detection.md）
- **快速近似搜索**：RandomProjectionIndexANN
- **小数据**：DummySNN 以获得精确结果

## 用例

### 模式匹配
查找一个模式在长序列中出现的位置：

<<<代码块_4>>>

### 主题发现
识别重复出现的模式：

<<<代码块_5>>>

### 时间序列检索
在数据库中查找相似的时间序列：

<<<代码块_6>>>

## 最佳实践

1. **窗口大小**：子序列方法的关键参数
   - 太小：捕获噪音
   - 太大：错过细粒度的图案
   - 经验法则：系列长度的 10-20%

2. **归一化**：大多数方法都假设 z 归一化子序列
   - 处理幅度变化
   - 注重形状相似性

3. **距离指标**：不同需求采用不同指标
   - 欧几里得：快速、基于形状
   - DTW：处理时间扭曲
   - 余弦：尺度不变

4. **排除区域**：对于主题发现，排除不重要的匹配
   - 通常设置为 0.5-1.0 × window_size
   - 防止发现重叠的情况

5. **性能**：
   - MASS 是 O(n log n) vs O(n²) 蛮力
   - 人工神经网络以准确性换取速度
   - GPU 加速可用于某些方法