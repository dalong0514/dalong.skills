<!-- 此文件由机器翻译自 best-practices.md -->

# Dask 最佳实践

## 性能优化原则

### 首先从更简单的解决方案开始

在使用 Dask 实施并行计算之前，请探索以下替代方案：
- 针对特定问题的更好算法
- 高效的文件格式（Parquet、HDF5、Zarr 而不是 CSV）
- 通过 Numba 或 Cython 编译代码
- 用于开发和测试的数据采样

这些替代方案通常比分布式系统提供更好的回报，并且应该在扩展到并行计算之前耗尽。

### 块大小策略

**关键规则**：块应该足够小，以便可以同时容纳许多块。

**推荐目标**：调整块大小，以便工作线程可以在每个核心容纳 10 个块，而不会超出可用内存。

**为什么重要**：
- 太大的块：内存溢出和低效的并行化
- 太小的块：过多的调度开销

**计算示例**：
- 8 核，32 GB RAM
- 目标：每个块约 400 MB（32 GB/8 核/10 块）

### 使用仪表板进行监控

Dask 仪表板提供了以下基本可见性：
- 工作状态和资源利用率
- 任务进度和瓶颈
- 内存使用模式
- 性能特点

访问仪表板以了解并行工作负载中实际缓慢的部分，而不是猜测优化情况。

## 要避免的关键陷阱

### 1. Dask之前不要在本地创建大对象

**错误的方法**：
```python
import pandas as pd
import dask.dataframe as dd

# Loads entire dataset into memory first
df = pd.read_csv('large_file.csv')
ddf = dd.from_pandas(df, npartitions=10)
```

**正确方法**：
<<<代码块_1>>>

**为什么**：使用 pandas 或 NumPy 加载数据首先会强制调度程序序列化这些对象并将其嵌入到任务图中，从而违背了并行计算的目的。

**关键原理**：使用Dask方法加载数据并使用Dask控制结果。

### 2. 避免重复调用compute()

**错误的方法**：
<<<代码块_2>>>

**正确方法**：
<<<代码块_3>>>

**为什么**：在循环中调用计算会阻止 Dask：
- 并行化不同的计算
- 共享中间结果
- 优化整体任务图

### 3. 不要构建过大的任务图

**症状**：
- 单次计算可处理数百万个任务
- 严重的调度开销
- 计算开始前的长时间延迟

**解决方案**：
- 增加块大小以减少任务数量
- 使用`map_partitions`或`map_blocks`来熔断操作
- 将计算分成较小的部分并保留中间内容
- 考虑问题是否真的需要分布式计算

**使用map_partitions的示例**：
<<<代码块_4>>>

## 基础设施注意事项

### 调度程序选择

**使用线程**：
- 使用 GIL 发布库（NumPy、Pandas、scikit-learn）进行数值工作
- 受益于共享内存的操作
- 具有数组/数据帧操作的单机工作负载

**使用流程**：
- 文本处理和Python集合操作
- GIL 绑定的纯 Python 代码
- 需要进程隔离的操作

**使用分布式调度程序**：
- 多机集群
- 需要诊断仪表板
- 异步API
- 更好的数据局部性处理

### 线程配置

**建议**：在数字工作负载上，每个进程的目标是大约 4 个线程。

**理由**：
- 并行性和开销之间的平衡
- 允许有效利用CPU核心
- 降低上下文切换成本

### 内存管理

**战略上坚持**：
<<<代码块_5>>>

**完成后清除内存**：
<<<代码块_6>>>

## 数据加载最佳实践

### 使用适当的文件格式

**对于表格数据**：
- Parquet：柱状、压缩、快速过滤
- CSV：仅适用于小数据或初始摄取

**对于数组数据**：
- HDF5：适用于数字数组
- Zarr：云原生、并行友好
- NetCDF：带有元数据的科学数据

### 优化数据摄取

**高效读取多个文件**：
```python
# Use glob patterns to read multiple files in parallel
ddf = dd.read_parquet('data/year=2024/month=*/day=*.parquet')
```

**尽早指定有用的列**：
```python
# Only read needed columns
ddf = dd.read_parquet('data.parquet', columns=['col1', 'col2', 'col3'])
```

## 常见模式和解决方案

### 模式：令人尴尬的并行问题

对于独立计算，请使用 Futures：
```python
from dask.distributed import Client

client = Client()
futures = [client.submit(func, arg) for arg in args]
results = client.gather(futures)
```

### 模式：数据预处理管道

使用 Bags 进行初始 ETL，然后转换为结构化格式：
```python
import dask.bag as db

# Process raw JSON
bag = db.read_text('logs/*.json').map(json.loads)
bag = bag.filter(lambda x: x['status'] == 'success')

# Convert to DataFrame for analysis
ddf = bag.to_dataframe()
```

### 模式：迭代算法

在迭代之间保留数据：
```python
data = dd.read_parquet('data.parquet')
data = data.persist()  # Keep in memory across iterations

for iteration in range(num_iterations):
    data = update_function(data)
    data = data.persist()  # Persist updated version
```

## 调试技巧

### 使用单线程调度程序

对于使用 pdb 进行调试或详细错误检查：
```python
import dask

dask.config.set(scheduler='synchronous')
result = computation.compute()  # Runs in single thread for debugging
```

### 检查任务图大小

计算之前，先检查一下任务数量：
```python
print(len(ddf.__dask_graph__()))  # Should be reasonable, not millions
```
### 首先验证小数据

在缩放之前测试小子集的逻辑：
```python
# Test on first partition
sample = ddf.head(1000)
# Validate results
# Then scale to full dataset
```

## 性能故障排除

### 症状：计算启动缓慢

**可能的原因**：任务图太大
**解决方案**：增加块大小或使用map_partitions

### 症状：内存错误

**可能的原因**：
- 块太大
- 中间结果太多
- 用户函数中的内存泄漏

**解决方案**：
- 减少块大小
- 策略性地使用 persist() 并在完成后删除
- 分析用户函数的内存问题

### 症状：并行化不佳

**可能的原因**：
- 数据依赖性阻碍并行性
- 块太大（没有足够的任务）
- GIL 与 Python 代码上的线程争用

**解决方案**：
- 重构计算以减少依赖性
- 增加分区数量
- 切换到 Python 代码的多处理调度程序