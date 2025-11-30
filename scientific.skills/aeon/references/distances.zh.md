<!-- 此文件由机器翻译自 distances.md -->

# 距离指标

Aeon 提供专门的距离函数来测量时间序列之间的相似性，与 aeon 和 scikit-learn 估计器兼容。

## 距离类别

### 弹性距离

允许系列之间灵活的时间对齐：

**动态时间扭曲系列：**
- `dtw` - 经典动态时间扭曲
- `ddtw` - 导数 DTW（比较导数）
- `wdtw` - 加权 DTW（按位置惩罚扭曲）
- `wddtw` - 加权导数 DTW
- `shape_dtw` - 基于形状的 DTW

**基于编辑：**
- `erp` - 使用真实惩罚编辑距离
- `edr` - 编辑真实序列上的距离
- `lcss` - 最长公共子序列
- `twe` - 时间扭曲编辑距离

**专业：**
- `msm` - 移动-拆分-合并距离
- `adtw` - 美国 DTW
- `sbd` - 基于形状的距离

**使用时间**：时间序列可能有时间变化、速度变化或相位差。

### 锁步距离

逐点比较时间序列而不对齐：

- `euclidean` - 欧氏距离（L2 范数）
- `manhattan` - 曼哈顿距离（L1 范数）
- `minkowski` - 广义 Minkowski 距离（Lp 范数）
- `squared` - 欧氏距离平方

**使用时间**：系列已对齐、需要计算速度或预计不会发生时间扭曲。

## 使用模式

### 计算单个距离

```python
from aeon.distances import dtw_distance

# Distance between two time series
distance = dtw_distance(x, y)

# With window constraint (Sakoe-Chiba band)
distance = dtw_distance(x, y, window=0.1)
```

### 成对距离矩阵

<<<代码块_1>>>

### 成本矩阵和对齐路径

<<<代码块_2>>>

### 与估算器一起使用

<<<代码块_3>>>

## 距离参数

### 窗口约束

限制翘曲路径偏差（提高速度并防止病态翘曲）：

<<<代码块_4>>>

### 标准化

控制在距离计算之前是否对序列进行 z 归一化：

<<<代码块_5>>>

### 距离特定参数

<<<代码块_6>>>

## 算法选择

### 按用例：

**时间错位**：DTW、DDTW、WDTW
**速度变化**：具有窗口约束的 DTW
**形状相似性**：形状 DTW、SBD
**编辑操作**：ERP、EDR、LCSS
**导数匹配**：DDTW
**计算速度**：欧几里得、曼哈顿
**异常稳健性**：曼哈顿、LCSS

### 按计算成本：

**最快**：欧几里得 (O(n))
**快速**：约束 DTW（O(nw)，其中 w 是窗口）
**中**：完整 DTW (O(n²))
**较慢**：复杂弹性距离（ERP、TWE、MSM）

## 快速参考表

|距离 |对齐|速度|稳健性|可解释性|
|----------|----------|--------|------------|--------------------|
|欧几里得|锁步 |非常快|低|高|
|大田 |弹性|中等|中等|中等|
| DDTW |弹性|中等|高|中等|
| WDTW |弹性|中等|中等|中等|
|企业资源规划|基于编辑 |慢|高|低|
| LCSS |基于编辑 |慢|非常高 |低|
|形状 DTW |弹性|中等|中等|高|

## 最佳实践

### 1. 标准化

大多数距离对比例敏感；适当时标准化：

```python
from aeon.transformations.collection import Normalizer

normalizer = Normalizer()
X_normalized = normalizer.fit_transform(X)
```

### 2. 窗口约束

对于 DTW 变体，使用窗口约束来提高速度和更好的泛化：

```python
# Start with 10-20% window
distance = dtw_distance(x, y, window=0.1)
```

### 3.系列长度

- 要求等长：大多数锁步距离
- 支持不等长：弹性距离（DTW、ERP等）

### 4.多元系列

大多数距离支持多元时间序列：

```python
# x.shape = (n_channels, n_timepoints)
distance = dtw_distance(x_multivariate, y_multivariate)
```

### 5.性能优化

- 使用 numba 编译的实现（aeon 中默认）
- 如果不需要对齐，请考虑锁步距离
- 使用窗口 DTW 而不是完整 DTW
- 预先计算距离矩阵以供重复使用

### 6. 选择合适的距离

```python
# Quick decision tree:
if series_aligned:
    use_distance = "euclidean"
elif need_speed:
    use_distance = "dtw"  # with window constraint
elif temporal_shifts_expected:
    use_distance = "dtw" or "shape_dtw"
elif outliers_present:
    use_distance = "lcss" or "manhattan"
elif derivatives_matter:
    use_distance = "ddtw" or "wddtw"
```

## 与 scikit-learn 集成

Aeon 距离与 sklearn 估计器一起使用：

```python
from sklearn.neighbors import KNeighborsClassifier
from aeon.distances import dtw_pairwise_distance

# Precompute distance matrix
X_train_distances = dtw_pairwise_distance(X_train)

# Use with sklearn
clf = KNeighborsClassifier(metric='precomputed')
clf.fit(X_train_distances, y_train)
```

## 可用的距离函数

获取所有可用距离的列表：

```python
from aeon.distances import get_distance_function_names

print(get_distance_function_names())
# ['dtw', 'ddtw', 'wdtw', 'euclidean', 'erp', 'edr', ...]
```

检索特定距离函数：

```python
from aeon.distances import get_distance_function

distance_func = get_distance_function("dtw")
result = distance_func(x, y, window=0.1)
```