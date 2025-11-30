<!-- 此文件由机器翻译自 performance.md -->

# 性能与优化

本参考资料涵盖了 Vaex 的性能特性，包括惰性求值、缓存、内存管理、异步操作以及处理海量数据集的优化策略。

## 理解惰性求值

惰性求值是 Vaex 性能的基础：

### 惰性求值是如何工作的

```python
import vaex

df = vaex.open('large_file.hdf5')

# No computation happens here - just defines what to compute
df['total'] = df.price * df.quantity
df['log_price'] = df.price.log()
mean_expr = df.total.mean()

# Computation happens here (when result is needed)
result = mean_expr  # Now the mean is actually calculated
```

**关键概念：**
- **表达式**是惰性的 - 它们定义计算而不执行它们
- **物化**在您访问结果时发生
- **查询优化**在执行前自动发生

### 评估何时进行？

<<<代码块_1>>>

## 批量操作延迟=True

一起执行多个操作以获得更好的性能：

### 基本延迟执行

<<<代码块_2>>>

### 多列延迟执行

<<<代码块_3>>>

### 何时使用delay=True

在以下情况下使用 `delay=True`：
- 计算多个聚合
- 计算多列的统计数据
- 构建仪表板或报告
- 任何需要多次传递数据的场景

<<<代码块_4>>>

## 异步操作

使用 async/await 异步处理数据：

### 异步与异步/等待

<<<代码块_5>>>

### 使用 Promise/Future

<<<代码块_6>>>

## 虚拟列与物化列

理解差异对于性能至关重要：

### 虚拟列（首选）

```python
# Virtual column - computed on-the-fly, zero memory
df['total'] = df.price * df.quantity
df['log_sales'] = df.sales.log()
df['full_name'] = df.first_name + ' ' + df.last_name

# Check if virtual
print(df.is_local('total'))  # False = virtual

# Benefits:
# - Zero memory overhead
# - Always up-to-date if source data changes
# - Fast to create
```

### 物化列

```python
# Materialize a virtual column
df['total_materialized'] = df['total'].values

# Or use materialize method
df = df.materialize(df['total'], inplace=True)

# Check if materialized
print(df.is_local('total_materialized'))  # True = materialized

# When to materialize:
# - Column computed repeatedly (amortize cost)
# - Complex expression used in many operations
# - Need to export data
```

### 决策：虚拟与实体化

```python
# Virtual is better when:
# - Column is simple (x + y, x * 2, etc.)
# - Column used infrequently
# - Memory is limited

# Materialize when:
# - Complex computation (multiple operations)
# - Used repeatedly in aggregations
# - Slows down other operations

# Example: Complex calculation used many times
df['complex'] = (df.x.log() * df.y.sqrt() + df.z ** 2).values  # Materialize
```

## 缓存策略

Vaex 会自动缓存一些操作，但您可以进一步优化：

### 自动缓存

```python
# First call computes and caches
mean1 = df.x.mean()  # Computes

# Second call uses cache
mean2 = df.x.mean()  # From cache (instant)

# Cache invalidated if DataFrame changes
df['new_col'] = df.x + 1
mean3 = df.x.mean()  # Recomputes
```

### 状态管理

```python
# Save DataFrame state (includes virtual columns)
df.state_write('state.json')

# Load state later
df_new = vaex.open('data.hdf5')
df_new.state_load('state.json')  # Restores virtual columns, selections
```

### 检查点模式

```python
# Export intermediate results for complex pipelines
df['processed'] = complex_calculation(df)

# Save checkpoint
df.export_hdf5('checkpoint.hdf5')

# Resume from checkpoint
df = vaex.open('checkpoint.hdf5')
# Continue processing...
```

## 内存管理

优化非常大的数据集的内存使用：

### 内存映射文件

```python
# HDF5 and Arrow are memory-mapped (optimal)
df = vaex.open('data.hdf5')  # No memory used until accessed

# File stays on disk, only accessed portions loaded to RAM
mean = df.x.mean()  # Streams through data, minimal memory
```

### 分块处理

```python
# Process large DataFrame in chunks
chunk_size = 1_000_000

for i1, i2, chunk in df.to_pandas_df(chunk_size=chunk_size):
    # Process chunk (careful: defeats Vaex's purpose)
    process_chunk(chunk)

# Better: Use Vaex operations directly (no chunking needed)
result = df.x.mean()  # Handles large data automatically
```

### 监控内存使用情况

```python
# Check DataFrame memory footprint
print(df.byte_size())  # Bytes used by materialized columns

# Check column memory
for col in df.get_column_names():
    if df.is_local(col):
        print(f"{col}: {df[col].nbytes / 1e9:.2f} GB")

# Profile operations
import vaex.profiler
with vaex.profiler():
    result = df.x.mean()
```

## 并行计算

Vaex 自动并行操作：

### 多线程

```python
# Vaex uses all CPU cores by default
import vaex

# Check/set thread count
print(vaex.multithreading.thread_count_default)
vaex.multithreading.thread_count_default = 8  # Use 8 threads

# Operations automatically parallelize
mean = df.x.mean()  # Uses all threads
```

### 使用 Dask 进行分布式计算

```python
# Convert to Dask for distributed processing
import vaex
import dask.dataframe as dd

# Create Vaex DataFrame
df_vaex = vaex.open('large_file.hdf5')

# Convert to Dask
df_dask = df_vaex.to_dask_dataframe()

# Process with Dask
result = df_dask.groupby('category')['value'].sum().compute()
```

## JIT 编译

Vaex 可以使用即时编译来进行自定义操作：

### 使用 Numba

```python
import vaex
import numba

# Define JIT-compiled function
@numba.jit
def custom_calculation(x, y):
    return x ** 2 + y ** 2

# Apply to DataFrame
df['custom'] = df.apply(custom_calculation,
                        arguments=[df.x, df.y],
                        vectorize=True)
```

### 自定义聚合

```python
@numba.jit
def custom_sum(a):
    total = 0
    for val in a:
        total += val * 2  # Custom logic
    return total

# Use in aggregation
result = df.x.custom_agg(custom_sum)
```

## 优化策略

### 策略 1：最小化物化

```python
# Bad: Creates many materialized columns
df['a'] = (df.x + df.y).values
df['b'] = (df.a * 2).values
df['c'] = (df.b + df.z).values

# Good: Keep virtual until final export
df['a'] = df.x + df.y
df['b'] = df.a * 2
df['c'] = df.b + df.z
# Only materialize if exporting:
# df.export_hdf5('output.hdf5')
```

### 策略 2：使用选择而不是过滤

```python
# Less efficient: Creates new DataFrames
df_high = df[df.value > 100]
df_low = df[df.value <= 100]
mean_high = df_high.value.mean()
mean_low = df_low.value.mean()

# More efficient: Use selections
df.select(df.value > 100, name='high')
df.select(df.value <= 100, name='low')
mean_high = df.value.mean(selection='high')
mean_low = df.value.mean(selection='low')
```

### 策略 3：批量聚合

```python
# Less efficient: Multiple passes
stats = {
    'mean': df.x.mean(),
    'std': df.x.std(),
    'min': df.x.min(),
    'max': df.x.max()
}

# More efficient: Single pass
delayed = [
    df.x.mean(delay=True),
    df.x.std(delay=True),
    df.x.min(delay=True),
    df.x.max(delay=True)
]
results = vaex.execute(delayed)
stats = dict(zip(['mean', 'std', 'min', 'max'], results))
```

### 策略 4：选择最佳文件格式

```python
# Slow: Large CSV
df = vaex.from_csv('huge.csv')  # Can take minutes

# Fast: HDF5 or Arrow
df = vaex.open('huge.hdf5')     # Instant
df = vaex.open('huge.arrow')    # Instant

# One-time conversion
df = vaex.from_csv('huge.csv', convert='huge.hdf5')
# Future loads: vaex.open('huge.hdf5')
```

### 策略 5：优化表达式

```python
# Less efficient: Repeated calculations
df['result'] = df.x.log() + df.x.log() * 2

# More efficient: Reuse calculations
df['log_x'] = df.x.log()
df['result'] = df.log_x + df.log_x * 2

# Even better: Combine operations
df['result'] = df.x.log() * 3  # Simplified math
```

## 性能分析

### 基本分析

```python
import time
import vaex

df = vaex.open('large_file.hdf5')

# Time operations
start = time.time()
result = df.x.mean()
elapsed = time.time() - start
print(f"Computed in {elapsed:.2f} seconds")
```

### 详细分析

```python
# Profile with context manager
with vaex.profiler():
    result = df.groupby('category').agg({'value': 'sum'})
# Prints detailed timing information
```

### 基准测试模式

```python
# Compare strategies
def benchmark_operation(operation, name):
    start = time.time()
    result = operation()
    elapsed = time.time() - start
    print(f"{name}: {elapsed:.3f}s")
    return result

# Test different approaches
benchmark_operation(lambda: df.x.mean(), "Direct mean")
benchmark_operation(lambda: df[df.x > 0].x.mean(), "Filtered mean")
benchmark_operation(lambda: df.x.mean(selection='positive'), "Selection mean")
```

## 常见性能问题及解决方案

### 问题：聚合缓慢

```python
# Problem: Multiple separate aggregations
for col in df.column_names:
    print(f"{col}: {df[col].mean()}")

# Solution: Batch with delay=True
delayed = [df[col].mean(delay=True) for col in df.column_names]
results = vaex.execute(delayed)
for col, result in zip(df.column_names, results):
    print(f"{col}: {result}")
```

### 问题：内存使用率高

```python
# Problem: Materializing large virtual columns
df['large_col'] = (complex_expression).values

# Solution: Keep virtual, or materialize and export
df['large_col'] = complex_expression  # Virtual
# Or: df.export_hdf5('with_new_col.hdf5')
```

### 问题：出口缓慢

```python
# Problem: Exporting with many virtual columns
df.export_csv('output.csv')  # Slow if many virtual columns

# Solution: Export to HDF5 or Arrow (faster)
df.export_hdf5('output.hdf5')
df.export_arrow('output.arrow')

# Or materialize first for CSV
df_materialized = df.materialize()
df_materialized.export_csv('output.csv')
```

### 问题：重复的复杂计算

```python
# Problem: Complex virtual column used repeatedly
df['complex'] = df.x.log() * df.y.sqrt() + df.z ** 3
result1 = df.groupby('cat1').agg({'complex': 'mean'})
result2 = df.groupby('cat2').agg({'complex': 'sum'})
result3 = df.complex.std()

# Solution: Materialize once
df['complex'] = (df.x.log() * df.y.sqrt() + df.z ** 3).values
# Or: df = df.materialize('complex')
```

## 性能最佳实践总结

1. **使用 HDF5 或 Arrow 格式** - 比 CSV 快几个数量级
2. **利用惰性计算** - 除非必要，否则不要强制计算
3. **延迟=True的批量操作** - 最小化数据传递
4. **保持列虚拟** - 仅在有益时才实现
5. **使用选择而不是过滤器** - 对于多个段更有效
6. **分析您的代码** - 在优化之前识别瓶颈
7. **避免 `.values` 和 `.to_pandas_df()`** - 在 Vaex 中保留操作
8. **自然并行化** - Vaex 自动使用所有内核
9. **导出为高效格式** - 检查点复杂管道
10. **优化表达式** - 简化数学并重用计算

## 相关资源

- 有关 DataFrame 基础知识：请参阅 `core_dataframes.md`
- 对于数据操作：参见`data_processing.md`
- 对于文件 I/O 优化：请参阅 `io_operations.md`