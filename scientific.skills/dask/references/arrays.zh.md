<!-- 此文件由机器翻译自 arrays.md -->

# Dask 数组

## 概述

Dask Array 使用阻塞算法实现 NumPy 的 ndarray 接口。它协调排列成网格的许多 NumPy 数组，以便利用跨多个核心的并行性，对大于可用内存的数据集进行计算。

## 核心理念

Dask 数组被分为块（块）：
- 每个块都是一个常规的 NumPy 数组
- 操作并行应用于每个块
- 结果自动合并
- 启用核外计算（数据大于 RAM）

## 关键能力

### Dask 阵列支持什么

**数学运算**：
- 算术运算（+、-、*、/）
- 标量函数（指数、对数、三角函数）
- 逐元素运算

**减少**：
- `sum()`、`mean()`、`std()`、`var()`
- 沿指定轴的减少
- `min()`、`max()`、`argmin()`、`argmax()`

**线性代数**：
- 张量收缩
- 点积和矩阵乘法
- 一些分解（SVD、QR）

**数据操作**：
- 换位
- 切片（标准和花式索引）
- 重塑
- 串联和堆叠

**阵列协议**：
- 通用函数（ufunc）
- 用于互操作性的 NumPy 协议

## 何时使用 Dask 数组

**在以下情况下使用 Dask 数组**：
- 数组超出可用 RAM
- 计算可以跨块并行化
- 使用 NumPy 风格的数值运算
- 需要将 NumPy 代码扩展到更大的数据集

**坚持使用 NumPy 时**：
- 数组可以轻松地放入内存中
- 运营需要全局数据视图
- 使用 Dask 中不可用的专用功能
- 单独使用 NumPy 性能就足够了

## 重要限制

Dask 数组故意不实现某些 NumPy 功能：

**未实施**：
- 大多数`np.linalg`函数（仅可用基本操作）
- 难以并行化的操作（如完全排序）
- 内存效率低下的操作（转换为列表，通过循环迭代）
- 许多专门功能（由社区需求驱动）

**解决方法**：对于不受支持的操作，请考虑将 `map_blocks` 与自定义 NumPy 代码结合使用。

## 创建 Dask 数组

### 来自 NumPy 数组
```python
import dask.array as da
import numpy as np

# Create from NumPy array with specified chunks
x = np.arange(10000)
dx = da.from_array(x, chunks=1000)  # Creates 10 chunks of 1000 elements each
```

### 随机数组
<<<代码块_1>>>

### 零、一和空
<<<代码块_2>>>

### 来自函数
<<<代码块_3>>>

### 从磁盘
<<<代码块_4>>>

## 常用操作

### 算术运算
<<<代码块_5>>>

### 减少
<<<代码块_6>>>

### 切片和索引
```python
# Standard slicing (returns Dask Array)
subset = x[1000:5000, 2000:8000]

# Fancy indexing
indices = [0, 5, 10, 15]
selected = x[indices, :]

# Boolean indexing
mask = x > 0.5
filtered = x[mask]
```

### 矩阵运算
```python
# Matrix multiplication
A = da.random.random((10000, 5000), chunks=(1000, 1000))
B = da.random.random((5000, 8000), chunks=(1000, 1000))
C = da.matmul(A, B)
result = C.compute()

# Dot product
dot_product = da.dot(A, B)

# Transpose
AT = A.T
```

### 线性代数
```python
# SVD (Singular Value Decomposition)
U, s, Vt = da.linalg.svd(A)
U_computed, s_computed, Vt_computed = dask.compute(U, s, Vt)

# QR decomposition
Q, R = da.linalg.qr(A)
Q_computed, R_computed = dask.compute(Q, R)

# Note: Only some linalg operations are available
```

### 重塑和操纵
```python
# Reshape
x = da.random.random((10000, 10000), chunks=(1000, 1000))
reshaped = x.reshape(5000, 20000)

# Transpose
transposed = x.T

# Concatenate
x1 = da.random.random((5000, 10000), chunks=(1000, 1000))
x2 = da.random.random((5000, 10000), chunks=(1000, 1000))
combined = da.concatenate([x1, x2], axis=0)

# Stack
stacked = da.stack([x1, x2], axis=0)
```

## 分块策略

分块对于 Dask Array 性能至关重要。

### 块大小指南

**良好的块尺寸**：
- 每个块：~10-100 MB（压缩）
- 每个数字数据块约 100 万个元素
- 并行性和开销之间的平衡

**计算示例**：
```python
# For float64 data (8 bytes per element)
# Target 100 MB chunks: 100 MB / 8 bytes = 12.5M elements

# For 2D array (10000, 10000):
x = da.random.random((10000, 10000), chunks=(1000, 1000))  # ~8 MB per chunk
```

### 查看块结构
```python
# Check chunks
print(x.chunks)  # ((1000, 1000, ...), (1000, 1000, ...))

# Number of chunks
print(x.npartitions)

# Chunk sizes in bytes
print(x.nbytes / x.npartitions)
```

### 重新分块
```python
# Change chunk sizes
x = da.random.random((10000, 10000), chunks=(500, 500))
x_rechunked = x.rechunk((2000, 2000))

# Rechunk specific dimension
x_rechunked = x.rechunk({0: 2000, 1: 'auto'})
```

## 使用 map_blocks 进行自定义操作

对于 Dask 中不可用的操作，请使用 `map_blocks`：

```python
import dask.array as da
import numpy as np

def custom_function(block):
    # Apply custom NumPy operation
    return np.fft.fft2(block)

x = da.random.random((10000, 10000), chunks=(1000, 1000))
result = da.map_blocks(custom_function, x, dtype=x.dtype)

# Compute
output = result.compute()
```

### 具有不同输出形状的map_blocks
```python
def reduction_function(block):
    # Returns scalar for each block
    return np.array([block.mean()])

result = da.map_blocks(
    reduction_function,
    x,
    dtype='float64',
    drop_axis=[0, 1],  # Output has no axes from input
    new_axis=0,        # Output has new axis
    chunks=(1,)        # One element per block
)
```

## 惰性评估和计算

### 惰性操作
```python
# All operations are lazy (instant, no computation)
x = da.random.random((10000, 10000), chunks=(1000, 1000))
y = x + 100
z = y.mean(axis=0)
result = z * 2

# Nothing computed yet, just task graph built
```

### 触发计算
```python
# Compute single result
final = result.compute()

# Compute multiple results efficiently
result1, result2 = dask.compute(operation1, operation2)
```

### 保留在内存中
```python
# Keep intermediate results in memory
x_cached = x.persist()

# Reuse cached results
y1 = (x_cached + 10).compute()
y2 = (x_cached * 2).compute()
```

## 保存结果

### 到 NumPy
```python
# Convert to NumPy (loads all in memory)
numpy_array = dask_array.compute()
```

### 到磁盘
```python
# Save to HDF5
import h5py
with h5py.File('output.hdf5', mode='w') as f:
    dset = f.create_dataset('/data', shape=x.shape, dtype=x.dtype)
    da.store(x, dset)

# Save to Zarr
import zarr
z = zarr.open('output.zarr', mode='w', shape=x.shape, dtype=x.dtype, chunks=x.chunks)
da.store(x, z)
```

## 性能考虑因素

### 高效运营
- 逐元素操作：非常高效
- 通过可并行操作减少：高效
- 沿块边界切片：高效
- 具有良好块对齐的矩阵运算：高效

### 昂贵的操作
- 跨多个块进行切片：需要数据移动
- 需要全局排序的操作：没有得到很好的支持
- 极其不规则的访问模式：性能差
- 块对齐较差的操作：需要重新分块

### 优化技巧

**1.选择合适的块尺寸**
```python
# Aim for balanced chunks
# Good: ~100 MB per chunk
x = da.random.random((100000, 10000), chunks=(10000, 10000))
```

**2.对齐操作块**
```python
# Make sure chunks align for operations
x = da.random.random((10000, 10000), chunks=(1000, 1000))
y = da.random.random((10000, 10000), chunks=(1000, 1000))  # Aligned
z = x + y  # Efficient
```

**3.使用适当的调度程序**
```python
# Arrays work well with threaded scheduler (default)
# Shared memory access is efficient
result = x.compute()  # Uses threads by default
```

**4.最大限度地减少数据传输**
```python
# Better: Compute on each chunk, then transfer results
means = x.mean(axis=1).compute()  # Transfers less data

# Worse: Transfer all data then compute
x_numpy = x.compute()
means = x_numpy.mean(axis=1)  # Transfers more data
```

## 常见模式

### 图像处理
```python
import dask.array as da

# Load large image stack
images = da.from_zarr('images.zarr')

# Apply filtering
def apply_gaussian(block):
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(block, sigma=2)

filtered = da.map_blocks(apply_gaussian, images, dtype=images.dtype)

# Compute statistics
mean_intensity = filtered.mean().compute()
```

### 科学计算
```python
# Large-scale numerical simulation
x = da.random.random((100000, 100000), chunks=(10000, 10000))

# Apply iterative computation
for i in range(num_iterations):
    x = da.exp(-x) * da.sin(x)
    x = x.persist()  # Keep in memory for next iteration

# Final result
result = x.compute()
```

### 数据分析
```python
# Load large dataset
data = da.from_zarr('measurements.zarr')

# Compute statistics
mean = data.mean(axis=0)
std = data.std(axis=0)
normalized = (data - mean) / std

# Save normalized data
da.to_zarr(normalized, 'normalized.zarr')
```

## 与其他工具集成

### XArray
```python
import xarray as xr
import dask.array as da

# XArray wraps Dask arrays with labeled dimensions
data = da.random.random((1000, 2000, 3000), chunks=(100, 200, 300))
dataset = xr.DataArray(
    data,
    dims=['time', 'y', 'x'],
    coords={'time': range(1000), 'y': range(2000), 'x': range(3000)}
)
```

### Scikit-learn（通过 Dask-ML）
```python
# Some scikit-learn compatible operations
from dask_ml.preprocessing import StandardScaler

X = da.random.random((10000, 100), chunks=(1000, 100))
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
```

## 调试技巧

### 可视化任务图
```python
# Visualize computation graph (for small arrays)
x = da.random.random((100, 100), chunks=(10, 10))
y = x + 1
y.visualize(filename='graph.png')
```

### 检查数组属性
```python
# Inspect before computing
print(f"Shape: {x.shape}")
print(f"Dtype: {x.dtype}")
print(f"Chunks: {x.chunks}")
print(f"Number of tasks: {len(x.__dask_graph__())}")
```

### 首先在小数组上进行测试
```python
# Test logic on small array
small_x = da.random.random((100, 100), chunks=(50, 50))
result_small = computation(small_x).compute()

# Validate, then scale
large_x = da.random.random((100000, 100000), chunks=(10000, 10000))
result_large = computation(large_x).compute()
```