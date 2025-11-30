<!-- 此文件由机器翻译自 distributed_computing.md -->

# 使用 Arboreto 进行分布式计算

Arboreto 利用 Dask 进行并行计算，实现从单机多核处理到多节点集群环境的高效 GRN 推理。

## 计算架构

GRN 推理本质上是可并行的：
- 每个目标基因的回归模型可以独立训练
- Arboreto 将计算表示为 Dask 任务图
- 任务分布在可用的计算资源上

## 本地多核处理（默认）

默认情况下，arboreto 使用本地计算机上所有可用的 CPU 核心：

```python
from arboreto.algo import grnboost2

# Automatically uses all local cores
network = grnboost2(expression_data=expression_matrix, tf_names=tf_names)
```

这对于大多数用例来说已经足够了，不需要额外的配置。

## 自定义本地 Dask 客户端

为了对本地资源进行细粒度控制，请创建自定义 Dask 客户端：

<<<代码块_1>>>

### 定制客户端的好处
- **资源控制**：限制CPU和内存使用
- **多次运行**：针对不同的参数集重复使用同一客户端
- **监控**：访问 Dask 仪表板以获取性能见解

## 使用同一客户端运行多个推理

重复使用单个 Dask 客户端进行具有不同参数的多次推理运行：

<<<代码块_2>>>

## 分布式集群计算

对于非常大的数据集，连接到在集群上运行的远程 Dask 分布式调度程序：

### 步骤 1：设置 Dask Scheduler（在簇头节点上）
<<<代码块_3>>>

### 步骤 2：启动 Dask Workers（在集群计算节点上）
<<<代码块_4>>>

### 第 3 步：从客户端连接
<<<代码块_5>>>

### 集群配置最佳实践

**工作人员配置**：
<<<代码块_6>>>

**对于大规模推理**：
- 使用更多具有中等内存的工作程序，而不是使用更少的具有大内存的工作程序
- 设置 `threads_per_worker=1` 以避免 scikit-learn 中的 GIL 争用
- 监控内存使用情况以防止worker被杀死

## 监控与调试

### Dask 仪表板

访问Dask仪表板进行实时监控：

```python
from distributed import Client

client = Client()  # Prints dashboard URL
# Dashboard available at: http://localhost:8787/status
```

仪表板显示：
- **任务进度**：已完成/待处理的任务数量
- **资源使用**：每个工作线程的 CPU、内存
- **任务流**：计算的实时可视化
- **性能**：瓶颈识别

### 详细输出

启用详细日志记录以跟踪推理进度：

```python
network = grnboost2(
    expression_data=expression_matrix,
    tf_names=tf_names,
    verbose=True
)
```

## 性能优化技巧

### 1. 数据格式
- **尽可能使用 Pandas DataFrame**：对于 Dask 操作比 NumPy 更高效
- **减少数据大小**：在推理之前过滤低方差基因

### 2.工作人员配置
- **CPU 密集型任务**：设置 `threads_per_worker=1`，增加 `n_workers`
- **内存绑定任务**：增加每个工作人员的 `memory_limit`

### 3. 集群设置
- **网络**：确保节点之间的高带宽、低延迟网络
- **存储**：对大型数据集使用共享文件系统或对象存储
- **调度**：分配专用节点以避免资源争用

### 4.转录因子过滤
- **限制 TF 列表**：提供特定的 TF 名称可减少计算量
```python
# Full search (slow)
network = grnboost2(expression_data=matrix)

# Filtered search (faster)
network = grnboost2(expression_data=matrix, tf_names=known_tfs)
```

## 示例：大规模单细胞分析

在集群上处理单细胞 RNA-seq 数据的完整工作流程：

```python
from distributed import Client
from arboreto.algo import grnboost2
import pandas as pd

if __name__ == '__main__':
    # Connect to cluster
    client = Client('tcp://cluster-scheduler:8786')

    # Load large single-cell dataset (50,000 cells x 20,000 genes)
    expression_data = pd.read_csv('scrnaseq_data.tsv', sep='\t')

    # Load cell-type-specific TFs
    tf_names = pd.read_csv('tf_list.txt', header=None)[0].tolist()

    # Run distributed inference
    network = grnboost2(
        expression_data=expression_data,
        tf_names=tf_names,
        client_or_address=client,
        verbose=True,
        seed=42
    )

    # Save results
    network.to_csv('grn_results.tsv', sep='\t', index=False)

    client.close()
```

这种方法可以分析在单台机器上不切实际的数据集。