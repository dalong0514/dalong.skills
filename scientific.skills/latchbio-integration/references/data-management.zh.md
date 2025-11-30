<!-- 此文件由机器翻译自 data-management.md -->

# 数据管理

## 概述
Latch 通过云存储抽象（LatchFile、LatchDir）和用于组织实验数据的结构化注册表系统提供全面的数据管理。

## 云存储：LatchFile 和 LatchDir

### 锁存文件

代表Latch云存储系统中的一个文件。

```python
from latch.types import LatchFile

# Create reference to existing file
input_file = LatchFile("latch:///data/sample.fastq")

# Access file properties
file_path = input_file.local_path  # Local path when executing
file_remote = input_file.remote_path  # Cloud storage path
```

### 锁存目录

代表Latch云存储系统中的一个目录。

<<<代码块_1>>>

### 路径格式

Latch 存储使用特殊的 URL 方案：
- **锁存路径**：`latch:///path/to/file`
- **本地路径**：工作流程执行期间自动解析
- **S3路径**：如果配置可以直接使用

### 文件传输

文件在本地执行和云存储之间自动传输：

<<<代码块_2>>>

### 全局模式

使用模式匹配查找文件：

<<<代码块_3>>>

## 注册系统

注册表提供了包含项目、表格和记录的结构化数据组织。

### 注册表架构

<<<代码块_4>>>

### 处理项目

<<<代码块_5>>>

### 使用表格

表存储结构化数据记录：

<<<代码块_6>>>

### 列类型

支持的数据类型：
- `string` - 文本数据
- `number` - 数值（整数或浮点数）
- `boolean` - 真/假值
- `date` - 日期值
- `file` - LatchFile 引用
- `directory` - LatchDir 引用
- `link` - 引用其他表中的记录
- `enum` - 预定义列表中的枚举值

### 使用记录

```python
from latch.registry.record import Record

# Create a record
record = Record.create(
    table_id=table.id,
    values={
        "sample_id": "S001",
        "condition": "treated",
        "replicate": 1,
        "fastq_file": LatchFile("latch:///data/S001.fastq")
    }
)

# Bulk create records
records = Record.bulk_create(
    table_id=table.id,
    records=[
        {"sample_id": "S001", "condition": "treated"},
        {"sample_id": "S002", "condition": "control"}
    ]
)

# Query records
all_records = Record.list(table_id=table.id)
filtered = Record.list(
    table_id=table.id,
    filter={"condition": "treated"}
)

# Update record
record.update(values={"replicate": 2})

# Delete record
record.delete()
```

### 链接记录

创建表之间的关系：

```python
# Define table with link column
results_table = Table.create(
    project_id=project.id,
    name="Results",
    columns=[
        {"name": "sample", "type": "link", "target_table": samples_table.id},
        {"name": "alignment_bam", "type": "file"},
        {"name": "gene_counts", "type": "file"}
    ]
)

# Create record with link
result_record = Record.create(
    table_id=results_table.id,
    values={
        "sample": sample_record.id,  # Link to sample record
        "alignment_bam": LatchFile("latch:///results/aligned.bam"),
        "gene_counts": LatchFile("latch:///results/counts.tsv")
    }
)

# Access linked data
sample_data = result_record.values["sample"].resolve()
```

### 枚举列

使用预定义值定义列：

```python
table = Table.create(
    project_id=project.id,
    name="Experiments",
    columns=[
        {
            "name": "status",
            "type": "enum",
            "options": ["pending", "running", "completed", "failed"]
        }
    ]
)
```

### 交易和批量更新

高效更新多条记录：

```python
from latch.registry.transaction import Transaction

# Start transaction
with Transaction() as txn:
    for record in records:
        record.update(values={"status": "processed"}, transaction=txn)
    # Changes committed when exiting context
```

## 与工作流程集成

### 在工作流程中使用注册表

```python
from latch import workflow, small_task
from latch.types import LatchFile
from latch.registry.table import Table
from latch.registry.record import Record

@small_task
def process_and_save(sample_id: str, table_id: str) -> str:
    # Get sample from registry
    table = Table.get(table_id=table_id)
    records = Record.list(
        table_id=table_id,
        filter={"sample_id": sample_id}
    )
    sample = records[0]

    # Process file
    input_file = sample.values["fastq_file"]
    # ... processing logic ...

    # Save results back to registry
    sample.update(values={
        "status": "completed",
        "results_file": output_file
    })

    return "Success"

@workflow
def registry_workflow(sample_id: str, table_id: str):
    """Workflow integrated with Registry"""
    return process_and_save(sample_id=sample_id, table_id=table_id)
```

### 数据自动执行工作流程

将工作流程配置为在将数据添加到注册表文件夹时自动运行：

```python
from latch.resources.launch_plan import LaunchPlan

# Create launch plan that watches a folder
launch_plan = LaunchPlan.create(
    workflow_name="rnaseq_pipeline",
    name="auto_process",
    trigger_folder="latch:///incoming_data",
    default_inputs={
        "output_dir": "latch:///results"
    }
)
```

## 帐户和工作空间管理

### 账户信息

```python
from latch.account import Account

# Get current account
account = Account.current()

# Account properties
workspace_id = account.id
workspace_name = account.name
```

### 团队工作空间

访问共享团队工作区：

```python
# List available workspaces
workspaces = Account.list()

# Switch workspace
Account.set_current(workspace_id="ws_789")
```

## 数据操作函数

### 连接数据

`latch.functions` 模块提供数据操作实用程序：

```python
from latch.functions import left_join, inner_join, outer_join, right_join

# Join tables
combined = left_join(
    left_table=table1,
    right_table=table2,
    on="sample_id"
)
```

### 过滤

```python
from latch.functions import filter_records

# Filter records
filtered = filter_records(
    table=table,
    condition=lambda record: record["replicate"] > 1
)
```

### 秘密管理

安全地存储和检索机密：

```python
from latch.functions import get_secret

# Retrieve secret in workflow
api_key = get_secret("api_key")
```

## 最佳实践

1. **路径组织**：使用一致的文件夹结构（例如，`/data`、`/results`、`/logs`）
2. **注册表架构**：在批量数据输入之前定义表架构
3. **链接记录**：使用链接来维护实验之间的关系
4. **批量操作**：使用事务更新多条记录
5. **文件命名**：使用一致的、描述性的文件命名约定
6. **元数据**：将实验元数据存储在注册表中以供追溯
7. **验证**：创建记录时验证数据类型
8. **清理**：定期归档或删除未使用的数据

## 常见模式

### 样本追踪

```python
# Create samples table
samples = Table.create(
    project_id=project.id,
    name="Samples",
    columns=[
        {"name": "sample_id", "type": "string"},
        {"name": "collection_date", "type": "date"},
        {"name": "raw_fastq_r1", "type": "file"},
        {"name": "raw_fastq_r2", "type": "file"},
        {"name": "status", "type": "enum", "options": ["pending", "processing", "complete"]}
    ]
)
```

### 结果组织

```python
# Create results table linked to samples
results = Table.create(
    project_id=project.id,
    name="Analysis Results",
    columns=[
        {"name": "sample", "type": "link", "target_table": samples.id},
        {"name": "alignment_bam", "type": "file"},
        {"name": "variants_vcf", "type": "file"},
        {"name": "qc_metrics", "type": "file"}
    ]
)
```