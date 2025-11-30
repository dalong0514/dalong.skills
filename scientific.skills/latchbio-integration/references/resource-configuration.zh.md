<!-- 此文件由机器翻译自 resource-configuration.md -->

# 资源配置

## 概述
Latch SDK 为工作流任务提供灵活的资源配置，从而能够在适当的计算基础设施（包括 CPU、GPU 和内存优化实例）上高效执行。

## 任务资源装饰器

### 标准装饰器

SDK 提供了预配置的任务装饰器来满足常见的资源需求：

#### @small_task
轻量级任务的默认配置：
```python
from latch import small_task

@small_task
def lightweight_processing():
    """Minimal resource requirements"""
    pass
```

**使用案例：**
- 文件解析和操作
- 简单的数据转换
- 快速质量控制检查
- 元数据操作

#### @large_task
增加 CPU 和内存以进行密集计算：
<<<代码块_1>>>

**使用案例：**
- 大文件处理
- 复杂的统计分析
- 组装任务
- 多线程操作

#### @small_gpu_task
以最少的资源启用 GPU：
<<<代码块_2>>>

**使用案例：**
- 神经网络推理
- 小规模机器学习预测
- GPU 加速库

#### @large_gpu_task
支持 GPU 并提供最大资源：
<<<代码块_3>>>

**使用案例：**
- 深度学习模型训练
- 蛋白质结构预测（AlphaFold）
- 大规模GPU计算

### 自定义任务配置

为了精确控制，请使用 `@custom_task` 装饰器：

<<<代码块_4>>>

#### 自定义任务参数

- **cpu**：CPU 核心数（整数）
- **内存**：内存（GB）（整数）
- **storage_gib**：GiB 中的临时存储（整数）
- **超时**：最大执行时间（以秒为单位）（整数）
- **gpu**：GPU 数量（整数，0 表示仅 CPU）
- **gpu_type**：特定 GPU 型号（字符串，例如“nvidia-tesla-v100”）

### 高级自定义配置

<<<代码块_5>>>

## GPU 配置

### GPU 类型

可用的 GPU 选项：
- **nvidia-tesla-k80**：用于测试的基本 GPU
- **nvidia-tesla-v100**：高性能训练
- **nvidia-tesla-a100**：最新一代的最高性能

### 多 GPU 任务

<<<代码块_6>>>

## 资源选择策略

### 按计算要求

**内存密集型任务：**
```python
@custom_task(cpu=4, memory=128)  # High memory, moderate CPU
def genome_assembly():
    pass
```

**CPU 密集型任务：**
```python
@custom_task(cpu=64, memory=32)  # High CPU, moderate memory
def parallel_alignment():
    pass
```

**I/O 密集型任务：**
```python
@custom_task(cpu=8, memory=16, storage_gib=1000)  # Large ephemeral storage
def data_preprocessing():
    pass
```

### 按工作流程阶段

**快速验证：**
```python
@small_task
def validate_inputs():
    """Fast input validation"""
    pass
```

**主要计算：**
```python
@large_task
def primary_analysis():
    """Resource-intensive analysis"""
    pass
```

**结果聚合：**
```python
@small_task
def aggregate_results():
    """Lightweight result compilation"""
    pass
```

## 工作流程资源规划

### 完整的管道示例

```python
from latch import workflow, small_task, large_task, large_gpu_task
from latch.types import LatchFile

@small_task
def quality_control(fastq: LatchFile) -> LatchFile:
    """QC doesn't need much resources"""
    return qc_output

@large_task
def alignment(fastq: LatchFile) -> LatchFile:
    """Alignment benefits from more CPU"""
    return bam_output

@large_gpu_task
def variant_calling(bam: LatchFile) -> LatchFile:
    """GPU-accelerated variant caller"""
    return vcf_output

@small_task
def generate_report(vcf: LatchFile) -> LatchFile:
    """Simple report generation"""
    return report

@workflow
def genomics_pipeline(input_fastq: LatchFile) -> LatchFile:
    """Resource-optimized genomics pipeline"""
    qc = quality_control(fastq=input_fastq)
    aligned = alignment(fastq=qc)
    variants = variant_calling(bam=aligned)
    return generate_report(vcf=variants)
```

## 超时配置

### 设置超时

```python
from latch import custom_task

@custom_task(
    cpu=8,
    memory=32,
    timeout=10800  # 3 hours in seconds
)
def long_running_analysis():
    """Analysis with extended timeout"""
    pass
```

### 超时最佳实践

1. **保守估计**：添加超出预期持续时间的缓冲时间
2. **监控实际运行时间**：根据实际执行数据进行调整
3. **默认超时**：大多数任务默认有1小时
4. **最大超时**：检查平台对超长作业的限制

## 存储配置

### 临时存储

配置中间文件的临时存储：

```python
@custom_task(
    cpu=8,
    memory=32,
    storage_gib=500  # 500 GB temporary storage
)
def process_large_dataset():
    """Task with large intermediate files"""
    # Ephemeral storage available at /tmp
    temp_file = "/tmp/intermediate_data.bam"
    pass
```

### 存储指南

- 默认存储通常足以满足大多数任务
- 为具有大型中间文件的任务指定更大的存储空间
- 任务完成后临时存储被清除
- 使用 LatchDir 满足持久存储需求

## 成本优化

### 资源效率技巧

1. **适当大小的资源**：不要过度分配
2. **使用适当的装饰器**：从标准装饰器开始
3. **仅在需要时使用GPU**：GPU任务成本更高
4. **并行小任务**：比一项大任务更好
5. **监控使用情况**：查看实际资源利用率

### 示例：经济高效的设计

```python
# INEFFICIENT: All tasks use large resources
@large_task
def validate_input():  # Over-provisioned
    pass

@large_task
def simple_transformation():  # Over-provisioned
    pass

# EFFICIENT: Right-sized resources
@small_task
def validate_input():  # Appropriate
    pass

@small_task
def simple_transformation():  # Appropriate
    pass

@large_task
def intensive_analysis():  # Appropriate
    pass
```

## 监控与调试

### 资源使用监控

在工作流程执行期间，监控：
- CPU利用率
- 内存使用情况
- GPU 利用率（如果适用）
- 执行持续时间
- 存储消耗

### 常见资源问题

**内存不足 (OOM)：**
```python
# Solution: Increase memory allocation
@custom_task(cpu=8, memory=64)  # Increased from 32 to 64 GB
def memory_intensive_task():
    pass
```

**超时：**
```python
# Solution: Increase timeout
@custom_task(cpu=8, memory=32, timeout=14400)  # 4 hours
def long_task():
    pass
```

**存储空间不足：**
```python
# Solution: Increase ephemeral storage
@custom_task(cpu=8, memory=32, storage_gib=1000)
def large_intermediate_files():
    pass
```

## 条件资源

根据输入动态分配资源：

```python
from latch import workflow, custom_task
from latch.types import LatchFile

def get_resource_config(file_size_gb: float):
    """Determine resources based on file size"""
    if file_size_gb < 10:
        return {"cpu": 4, "memory": 16}
    elif file_size_gb < 100:
        return {"cpu": 16, "memory": 64}
    else:
        return {"cpu": 32, "memory": 128}

# Note: Resource decorators must be static
# Use multiple task variants for different sizes
@custom_task(cpu=4, memory=16)
def process_small(file: LatchFile) -> LatchFile:
    pass

@custom_task(cpu=16, memory=64)
def process_medium(file: LatchFile) -> LatchFile:
    pass

@custom_task(cpu=32, memory=128)
def process_large(file: LatchFile) -> LatchFile:
    pass
```

## 最佳实践总结

1. **从小处开始**：从标准装饰器开始，根据需要扩大规模
2. **先分析**：运行测试执行以确定实际需求
3. **谨慎使用 GPU**：仅在算法支持时才使用 GPU
4. **并行设计**：尽可能分解成更小的任务
5. **监控和调整**：审查执行指标并优化
6. **文件要求**：评论为什么需要特定资源
7. **本地测试**：注册前在本地使用Docker进行验证
8. **考虑成本**：平衡性能与成本效率

## 特定于平台的注意事项

### 可用资源

Latch 平台提供：
- CPU 实例：最多 96 个核心
- 内存：高达 768 GB
- GPU：K80、V100、A100 变体
- 存储：可配置的临时存储

### 资源限制

检查当前平台限制：
- 每个任务的最大CPU数
- 每个任务的最大内存
- 最大GPU分配
- 最大并发任务数

### 配额和限制

请注意工作空间配额：
- 并发执行总数
- 总资源分配
- 存储限制
- 执行时间限制

如果需要，请联系 Latch 支持人员以增加配额。