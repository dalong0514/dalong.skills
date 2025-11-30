<!-- 此文件由机器翻译自 api_reference.md -->

# FlowIO API 参考

## 概述

FlowIO 是一个用于读写流式细胞术标准 (FCS) 文件的 Python 库。它以最小的依赖性支持 FCS 版本 2.0、3.0 和 3.1。

## 安装

```bash
pip install flowio
```

支持Python 3.9及更高版本。

## 核心课程

### 流量数据

用于处理 FCS 文件的主要类。

#### 构造函数

<<<代码块_1>>>

**参数：**
- `fcs_file`：文件路径 (str)、路径对象或文件句柄
- `ignore_offset_error` (bool)：忽略偏移错误（默认值：False）
- `ignore_offset_discrepancy` (bool)：忽略 HEADER 和 TEXT 部分之间的偏移量差异（默认值：False）
- `use_header_offsets` (bool)：使用 HEADER 部分偏移量而不是 TEXT 部分（默认值：False）
- `only_text` (bool): 只解析 TEXT 段，跳过 DATA 和 ANALYSIS (默认值: False)
- `nextdata_offset` (int): 读取多数据集文件的字节偏移量
- `null_channel_list`（列表）：要排除的空通道的 PnN 标签列表

#### 属性

**文件信息：**
- `name`：FCS 文件的名称
- `file_size`：文件大小（以字节为单位）
- `version`：FCS 版本（例如“3.0”、“3.1”）
- `header`：包含 HEADER 段信息的字典
- `data_type`：数据格式类型（'I'、'F'、'D'、'A'）

**频道信息：**
- `channel_count`：数据集中的通道数
- `channels`：将通道编号映射到通道信息的字典
- `pnn_labels`：PnN（短通道名称）标签列表
- `pns_labels`：PnS（描述性染色名称）标签列表
- `pnr_values`：每个通道的 PnR（范围）值列表
- `fluoro_indices`：荧光通道索引列表
- `scatter_indices`：分散通道的索引列表
- `time_index`：时间通道的索引（或无）
- `null_channels`：空通道索引列表

**事件数据：**
- `event_count`：数据集中的事件数（行）
- `events`：原始事件数据（字节）

**元数据：**
- `text`：TEXT 段键值对字典
- `analysis`：ANALYSIS 段键值对字典（如果存在）

#### 方法

##### as_array()

<<<代码块_2>>>

以二维 NumPy 数组形式返回事件数据。

**参数：**
- `preprocess` (bool)：应用增益、对数和时间缩放变换（默认值：True）

**退货：**
- NumPy ndarray 形状（event_count，channel_count）

**示例：**
<<<代码块_3>>>

##### write_fcs()

<<<代码块_4>>>

将 FlowData 实例导出为新的 FCS 文件。

**参数：**
- `filename` (str): 输出文件路径
- `metadata` (dict)：要添加/更新的文本段关键字的可选字典

**示例：**
<<<代码块_5>>>

**注意：** 导出为带有单精度浮点数据的 FCS 3.1。

## 实用函数

### read_multiple_data_sets()

<<<代码块_6>>>

从包含多个数据集的 FCS 文件中读取所有数据集。

**参数：**
- 与 FlowData 构造函数相同（`nextdata_offset` 除外）

**退货：**
- FlowData 实例列表，每个数据集一个

**示例：**
```python
from flowio import read_multiple_data_sets

datasets = read_multiple_data_sets('multi_dataset.fcs')
print(f"Found {len(datasets)} datasets")
for i, dataset in enumerate(datasets):
    print(f"Dataset {i}: {dataset.event_count} events")
```

### create_fcs()

```python
create_fcs(filename,
           event_data,
           channel_names,
           opt_channel_names=None,
           metadata=None)
```

从事件数据创建新的 FCS 文件。

**参数：**
- `filename` (str): 输出文件路径
- `event_data` (ndarray)：事件数据的二维 NumPy 数组（行=事件，列=通道）
- `channel_names`（列表）：PnN（短）通道名称列表
- `opt_channel_names`（列表）：PnS（描述性）通道名称的可选列表
- `metadata` (dict)：TEXT 段关键字的可选字典

**示例：**
```python
import numpy as np
from flowio import create_fcs

# Create synthetic data
events = np.random.rand(10000, 5)
channels = ['FSC-A', 'SSC-A', 'FL1-A', 'FL2-A', 'Time']
opt_channels = ['Forward Scatter', 'Side Scatter', 'FITC', 'PE', 'Time']

create_fcs('synthetic.fcs',
           events,
           channels,
           opt_channel_names=opt_channels,
           metadata={'$SRC': 'Synthetic data'})
```

## 异常类

### FlowIO警告

非关键问题的通用警告类别。

### PnE警告

当 FCS 文件创建期间 PnE 值无效时引发警告。

### FlowIOException

FlowIO 错误的基本异常类。

### FCS解析错误

当解析 FCS 文件出现问题时引发。

### 数据偏移差异错误

当 HEADER 和 TEXT 部分为数据段提供不同的字节偏移量时引发。

**解决方法：** 创建 FlowData 实例时使用 `ignore_offset_discrepancy=True` 参数。

### 多个数据集错误

尝试使用标准 FlowData 构造函数读取包含多个数据集的文件时引发。

**解决方案：** 使用 `read_multiple_data_sets()` 函数代替。

## FCS 文件结构参考

FCS 文件由四个部分组成：
1. **HEADER**：包含FCS版本和其他段的字节位置
2. **TEXT**：键值元数据对（分隔格式）
3. **DATA**：原始事件数据（二进制、浮点或 ASCII）
4. **分析**（可选）：数据处理结果

### 常见文本段关键字

- `$BEGINDATA`、`$ENDDATA`：DATA 段的字节偏移量
- `$BEGINANALYSIS`、`$ENDANALYSIS`：ANALYSIS 段的字节偏移量
- `$BYTEORD`：字节顺序（1,2,3,4 表示小端；4,3,2,1 表示大端）
- `$DATATYPE`：数据类型（'I'=整数，'F'=浮点，'D'=双精度，'A'=ASCII）
- `$MODE`：数据模式（'L'=列表模式，最常见）
- `$NEXTDATA`：到下一个数据集的偏移量（如果是单个数据集则为 0）
- `$PAR`：参数数量（通道）
- `$TOT`：事件总数
- `PnN`：参数 n 的简称
- `PnS`：参数 n 的描述性污点名称
- `PnR`：参数 n 的范围（最大值）
- `PnE`：参数 n 的放大指数（格式：“a,b”，其中值 = a * 10^(b*x)）
- `PnG`：参数 n 的放大增益

## 通道类型

FlowIO 自动对通道进行分类：

- **散射通道**：FSC（前向散射）、SSC（侧向散射）
- **荧光通道**：FL1、FL2、FITC、PE 等。
- **时间通道**：通常标记为“时间”

通过以下方式访问索引：
- `flow_data.scatter_indices`
- `flow_data.fluoro_indices`
- `flow_data.time_index`

## 数据预处理

当调用 `as_array(preprocess=True)` 时，FlowIO 适用：

1. **增益缩放**：乘以 PnG 值
2. **对数变换**：应用 PnE 指数变换（如果存在）
3. **时间缩放**：将时间值转换为适当的单位

要访问原始的、未处理的数据：`as_array(preprocess=False)`

## 最佳实践

1. **内存效率**：仅需要元数据时使用`only_text=True`
2. **错误处理**：将文件操作包装在 FCSParsingError 的 try- except 块中
3. **多数据集文件**：如果不确定数据集计数，请始终使用 `read_multiple_data_sets()`
4. **偏移问题**：如果遇到偏移错误，请尝试`ignore_offset_discrepancy=True`
5. **通道选择**：解析时使用null_channel_list排除不需要的通道

## 与 FlowKit 集成

对于高级流式细胞术分析（包括补偿、门控和 GatingML 支持），请考虑将 FlowKit 库与 FlowIO 一起使用。 FlowKit 提供了构建在 FlowIO 文件解析功能之上的更高级别的抽象。

## 工作流程示例

### 基本文件读取

```python
from flowio import FlowData

# Read FCS file
flow = FlowData('experiment.fcs')

# Print basic info
print(f"Version: {flow.version}")
print(f"Events: {flow.event_count}")
print(f"Channels: {flow.channel_count}")
print(f"Channel names: {flow.pnn_labels}")

# Get event data
events = flow.as_array()
print(f"Data shape: {events.shape}")
```

### 元数据提取

```python
from flowio import FlowData

flow = FlowData('sample.fcs', only_text=True)

# Access metadata
print(f"Acquisition date: {flow.text.get('$DATE', 'N/A')}")
print(f"Instrument: {flow.text.get('$CYT', 'N/A')}")

# Channel information
for i, (pnn, pns) in enumerate(zip(flow.pnn_labels, flow.pns_labels)):
    print(f"Channel {i}: {pnn} ({pns})")
```

### 创建新的 FCS 文件

```python
import numpy as np
from flowio import create_fcs

# Generate or process data
data = np.random.rand(5000, 3) * 1000

# Define channels
channels = ['FSC-A', 'SSC-A', 'FL1-A']
stains = ['Forward Scatter', 'Side Scatter', 'GFP']

# Create FCS file
create_fcs('output.fcs',
           data,
           channels,
           opt_channel_names=stains,
           metadata={
               '$SRC': 'Python script',
               '$DATE': '19-OCT-2025'
           })
```

### 处理多数据集文件

```python
from flowio import read_multiple_data_sets

# Read all datasets
datasets = read_multiple_data_sets('multi.fcs')

# Process each dataset
for i, dataset in enumerate(datasets):
    print(f"\nDataset {i}:")
    print(f"  Events: {dataset.event_count}")
    print(f"  Channels: {dataset.pnn_labels}")

    # Get data array
    events = dataset.as_array()
    mean_values = events.mean(axis=0)
    print(f"  Mean values: {mean_values}")
```

### 修改并重新导出

```python
from flowio import FlowData

# Read original file
flow = FlowData('original.fcs')

# Get event data
events = flow.as_array(preprocess=False)

# Modify data (example: apply custom transformation)
events[:, 0] = events[:, 0] * 1.5  # Scale first channel

# Note: Currently, FlowIO doesn't support direct modification of event data
# For modifications, use create_fcs() instead:
from flowio import create_fcs

create_fcs('modified.fcs',
           events,
           flow.pnn_labels,
           opt_channel_names=flow.pns_labels,
           metadata=flow.text)
```