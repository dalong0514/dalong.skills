<!-- 此文件由机器翻译自 python-api.md -->

# Python API 参考

gtars Python 绑定的综合参考。

## 安装

```bash
# Install gtars Python package
uv pip install gtars

# Or with pip
pip install gtars
```

## 核心课程

### 区域集

管理基因组区间的集合：

<<<代码块_1>>>

### 区域运营

对区域集执行操作：

<<<代码块_2>>>

### 设置操作

对基因组区域执行集合操作：

<<<代码块_3>>>

## 数据导出

### 写入 BED 文件

将区域导出为 BED 格式：

<<<代码块_4>>>

### 格式转换

格式之间的转换：

<<<代码块_5>>>

## NumPy 集成

与 NumPy 数组无缝集成：

<<<代码块_6>>>

## 并行处理

利用大型数据集的并行处理：

```python
# Enable parallel processing
regions = gtars.RegionSet.from_bed("large_file.bed", parallel=True)

# Parallel operations
result = regions.parallel_apply(custom_function)
```

## 内存管理

大型数据集的高效内存使用：

```python
# Stream large BED files
for chunk in gtars.RegionSet.stream_bed("large_file.bed", chunk_size=10000):
    process_chunk(chunk)

# Memory-mapped mode
regions = gtars.RegionSet.from_bed("large_file.bed", mmap=True)
```

## 错误处理

处理常见错误：

```python
try:
    regions = gtars.RegionSet.from_bed("file.bed")
except gtars.FileNotFoundError:
    print("File not found")
except gtars.InvalidFormatError as e:
    print(f"Invalid BED format: {e}")
except gtars.ParseError as e:
    print(f"Parse error at line {e.line}: {e.message}")
```

## 配置

配置 gtar 行为：

```python
# Set global options
gtars.set_option("parallel.threads", 4)
gtars.set_option("memory.limit", "4GB")
gtars.set_option("warnings.strict", True)

# Context manager for temporary options
with gtars.option_context("parallel.threads", 8):
    # Use 8 threads for this block
    regions = gtars.RegionSet.from_bed("large_file.bed", parallel=True)
```

## 日志记录

启用日志记录以进行调试：

```python
import logging

# Enable gtars logging
gtars.set_log_level("DEBUG")

# Or use Python logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("gtars")
```

## 性能提示

- 对大型数据集使用并行处理
- 为非常大的文件启用内存映射模式
- 尽可能流式传输数据以减少内存使用
- 适用时，在操作前对区域进行预排序
- 使用 NumPy 数组进行数值计算
- 缓存经常访问的数据