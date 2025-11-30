<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：达斯克
描述：“并行/分布式计算。将 pandas/NumPy 扩展到内存之外、并行数据帧/数组、多文件处理、任务图，适用于大于 RAM 的数据集和并行工作流程。”
---

# 达斯克

## 概述

Dask 是一个用于并行和分布式计算的 Python 库，它支持三个关键功能：
- 对于超过可用 RAM 的数据，在单台机器上执行**大于内存的执行**
- **并行处理**可提高跨多个核心的计算速度
- **分布式计算**支持跨多台机器的 TB 级数据集

Dask 可从笔记本电脑（处理约 100 GiB）扩展到集群（处理约 100 TiB），同时保持熟悉的 Python API。

## 何时使用此技能

该技能应该在以下情况下使用：
- 处理超出可用 RAM 的数据集
- 将 pandas 或 NumPy 操作扩展到更大的数据集
- 并行计算以提高性能
- 高效处理多个文件（CSV、Parquet、JSON、文本日志）
- 构建具有任务依赖性的自定义并行工作流程
- 在多个核心或机器之间分配工作负载

## 核心能力

Dask 提供五个主要组件，每个组件适合不同的用例：

### 1. DataFrames - 并行 Pandas 操作

**目的**：通过并行处理将 pandas 操作扩展到更大的数据集。

**何时使用**：
- 表格数据超出可用 RAM
- 需要一起处理多个 CSV/Parquet 文件
- Pandas 操作速度慢，需要并行化
- 从 pandas 原型扩展到生产

**参考文档**：有关 Dask DataFrames 的全面指南，请参阅 `references/dataframes.md`，其中包括：
- 读取数据（单个文件、多个文件、glob 模式）
- 常用操作（过滤、分组、连接、聚合）
- 使用 `map_partitions` 进行自定义操作
- 性能优化技巧
- 常见模式（ETL、时间序列、多文件处理）

**简单示例**：
```python
import dask.dataframe as dd

# Read multiple files as single DataFrame
ddf = dd.read_csv('data/2024-*.csv')

# Operations are lazy until compute()
filtered = ddf[ddf['value'] > 100]
result = filtered.groupby('category').mean().compute()
```

**要点**：
- 操作是惰性的（构建任务图），直到`.compute()`被调用
- 使用`map_partitions`进行高效的自定义操作
- 在处理来自其他来源的结构化数据时尽早转换为 DataFrame

### 2. 数组 - 并行 NumPy 运算

**目的**：使用阻塞算法将 NumPy 功能扩展到大于内存的数据集。

**何时使用**：
- 数组超出可用 RAM
- NumPy 操作需要并行化
- 使用科学数据集（HDF5、Zarr、NetCDF）
- 需要并行线性代数或数组运算

**参考文档**：有关 Dask 数组的综合指南，请参阅 `references/arrays.md`，其中包括：
- 创建数组（从 NumPy、随机、从磁盘）
- 分块策略和优化
- 常见运算（算术、归约、线性代数）
- 使用 `map_blocks` 进行自定义操作
- 与 HDF5、Zarr 和 XArray 集成

**简单示例**：
<<<代码块_1>>>

**要点**：
- 块大小至关重要（目标是每个块约 100 MB）
- 对块的操作并行进行
- 需要时重新分块数据以实现高效运营
- 使用 `map_blocks` 进行 Dask 中不可用的操作

### 3. Bags - 非结构化数据的并行处理

**目的**：通过函数操作处理非结构化或半结构化数据（文本、JSON、日志）。

**何时使用**：
- 处理文本文件、日志或 JSON 记录
- 结构化分析前的数据清理和ETL
- 使用不适合数组/数据帧格式的Python对象
- 需要内存高效的流处理

**参考文档**：有关 Dask Bags 的综合指南，请参阅 `references/bags.md`，其中包括：
- 读取文本和 JSON 文件
- 功能操作（映射、过滤、折叠、分组）
- 转换为数据帧
- 常见模式（日志分析、JSON处理、文本处理）
- 性能考虑

**简单示例**：
<<<代码块_2>>>

**要点**：
- 用于初始数据清理，然后转换为DataFrame/Array
- 使用 `foldby` 而不是 `groupby` 以获得更好的性能
- 操作是流式的并且内存效率高
- 转换为结构化格式（DataFrame）以进行复杂操作

### 4. Futures - 基于任务的并行化

**目的**：通过对任务执行和依赖关系的细粒度控制来构建自定义并行工作流程。

**何时使用**：
- 构建动态、不断发展的工作流程
- 需要立即执行任务（不是懒惰）
- 计算取决于运行时条件
- 实现自定义并行算法
- 需要有状态计算

**参考文档**：有关 Dask Futures 的综合指南，请参阅 `references/futures.md`，其中包括：
- 设置分布式客户端
- 提交任务并与 future 合作
- 任务依赖性和数据移动
- 高级协调（队列、锁、事件、参与者）
- 常见模式（参数扫描、动态任务、迭代算法）

**简单示例**：
<<<代码块_3>>>

**要点**：
- 需要分布式客户端（即使是单机）
- 任务提交后立即执行
- 大数据预分散，避免重复传输
- 每个任务约 1 毫秒的开销（不适合数百万个小任务）
- 使用参与者进行有状态的工作流程

### 5. 调度程序 - 执行后端

**用途**：控制 Dask 任务的执行方式和位置（线程、进程、分布式）。

**何时选择调度程序**：
- **线程**（默认）：NumPy/Pandas 操作、GIL 释放库、共享内存优势
- **进程**：纯Python代码、文本处理、GIL绑定操作
- **同步**：使用 pdb 进行调试、分析、理解错误
- **分布式**：需要仪表板、多机集群、高级功能

**参考文档**：有关 Dask Scheduler 的全面指南，请参阅 `references/schedulers.md`，其中包括：
- 详细的调度程序描述和特性
- 配置方法（全局、上下文管理器、每次计算）
- 性能考虑和开销
- 常见模式和故障排除
- 线程配置以获得最佳性能

**简单示例**：
<<<代码块_4>>>

**要点**：
- 线程：最低开销（~10 µs/任务），最适合数字工作
- 进程：避免 GIL（~10 毫秒/任务），最适合 Python 工作
- 分布式：监控仪表板（~1 毫秒/任务），可扩展到集群
- 可以根据计算或全局切换调度程序

## 最佳实践

有关全面的性能优化指南、内存管理策略以及要避免的常见陷阱，请参阅`references/best-practices.md`。主要原则包括：

### 从更简单的解决方案开始
在使用 Dask 之前，请探索：
- 更好的算法
- 高效的文件格式（Parquet 而不是 CSV）
- 编译代码（Numba、Cython）
- 数据采样

### 关键性能规则

**1.不要在本地加载数据然后交给 Dask**
<<<代码块_5>>>

**2.避免重复的compute()调用**
<<<代码块_6>>>

**3.不要构建过大的任务图**
- 如果有数百万个任务，则增加块大小
- 使用`map_partitions`/`map_blocks`来熔断操作
- 检查任务图大小：`len(ddf.__dask_graph__())`

**4.选择合适的块大小**
- 目标：每个块约 100 MB（或工作内存中每个核心 10 个块）
- 太大：内存溢出
- 太小：调度开销

**5.使用仪表板**
```python
from dask.distributed import Client
client = Client()
print(client.dashboard_link)  # Monitor performance, identify bottlenecks
```

## 常见工作流程模式

### ETL 管道
```python
import dask.dataframe as dd

# Extract: Read data
ddf = dd.read_csv('raw_data/*.csv')

# Transform: Clean and process
ddf = ddf[ddf['status'] == 'valid']
ddf['amount'] = ddf['amount'].astype('float64')
ddf = ddf.dropna(subset=['important_col'])

# Load: Aggregate and save
summary = ddf.groupby('category').agg({'amount': ['sum', 'mean']})
summary.to_parquet('output/summary.parquet')
```

### 非结构化到结构化管道
```python
import dask.bag as db
import json

# Start with Bag for unstructured data
bag = db.read_text('logs/*.json').map(json.loads)
bag = bag.filter(lambda x: x['status'] == 'valid')

# Convert to DataFrame for structured analysis
ddf = bag.to_dataframe()
result = ddf.groupby('category').mean().compute()
```

### 大规模数组计算
```python
import dask.array as da

# Load or create large array
x = da.from_zarr('large_dataset.zarr')

# Process in chunks
normalized = (x - x.mean()) / x.std()

# Save result
da.to_zarr(normalized, 'normalized.zarr')
```

### 自定义并行工作流程
```python
from dask.distributed import Client

client = Client()

# Scatter large dataset once
data = client.scatter(large_dataset)

# Process in parallel with dependencies
futures = []
for param in parameters:
    future = client.submit(process, data, param)
    futures.append(future)

# Gather results
results = client.gather(futures)
```

## 选择正确的组件

使用此决策指南来选择合适的 Dask 组件：

**数据类型**：
- 表格数据 → **数据帧**
- 数字数组 → **数组**
- 文本/JSON/日志 → **Bags** （然后转换为 DataFrame）
- 自定义 Python 对象 → **Bags** 或 **Futures**

**操作类型**：
- 标准 pandas 操作 → **DataFrames**
- 标准 NumPy 运算 → **数组**
- 自定义并行任务 → **Futures**
- 文本处理/ETL → **包**

**控制级别**：
- 高级、自动 → **数据帧/数组**
- 低级，手动 → **期货**

**工作流程类型**：
- 静态计算图 → **DataFrames/Arrays/Bags**
- 动态、不断发展 → **未来**

## 集成注意事项

### 文件格式
- **高效**：Parquet、HDF5、Zarr（柱状、压缩、并行友好）
- **兼容但速度较慢**：CSV（仅用于初始摄取）
- **对于阵列**：HDF5、Zarr、NetCDF

### 集合之间的转换
```python
# Bag → DataFrame
ddf = bag.to_dataframe()

# DataFrame → Array (for numeric data)
arr = ddf.to_dask_array(lengths=True)

# Array → DataFrame
ddf = dd.from_dask_array(arr, columns=['col1', 'col2'])
```

### 与其他库
- **XArray**：用标记尺寸（地理空间、成像）包装 Dask 数组
- **Dask-ML**：使用 scikit-learn 兼容 API 进行机器学习
- **分布式**：高级集群管理和监控

## 调试与开发

### 迭代开发工作流程

1. **使用同步调度程序对小数据进行测试**：
```python
dask.config.set(scheduler='synchronous')
result = computation.compute()  # Can use pdb, easy debugging
```

2. **使用示例上的线程进行验证**：
```python
sample = ddf.head(1000)  # Small sample
# Test logic, then scale to full dataset
```
3. **分布式监控规模**：
```python
from dask.distributed import Client
client = Client()
print(client.dashboard_link)  # Monitor performance
result = computation.compute()
```

### 常见问题

**内存错误**：
- 减少块大小
- 策略性地使用`persist()`并在完成后删除
- 检查自定义函数中的内存泄漏

**慢启动**：
- 任务图太大（增加块大小）
- 使用`map_partitions`或`map_blocks`来减少任务

**并行化较差**：
- 块太大（增加分区数量）
- 在 Python 代码中使用线程（切换到进程）
- 数据依赖性阻碍并行性

## 参考文件

可以根据需要阅读所有参考文档文件以获取详细信息：

- `references/dataframes.md` - 完整的 Dask DataFrame 指南
- `references/arrays.md` - 完整的 Dask 阵列指南
- `references/bags.md` - 完整 Dask Bag 指南
- `references/futures.md` - 完整的 Dask Futures 和分布式计算指南
- `references/schedulers.md` - 完整的调度程序选择和配置指南
- `references/best-practices.md` - 全面的性能优化和故障排除

当用户需要有关特定 Dask 组件、操作或模式的详细信息（超出此处提供的快速指南）时，请加载这些文件。