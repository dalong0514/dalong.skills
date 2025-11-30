<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 基准测试集成
描述：“Benchling 研发平台集成。通过 API 访问注册表（DNA、蛋白质）、库存、ELN 条目、工作流程、构建 Benchling 应用程序、查询数据仓库，实现实验室数据管理自动化。”
---

# 基准测试集成

## 概述

Benchling 是一个生命科学研发云平台。通过 Python SDK 和 REST API 以编程方式访问注册表实体（DNA、蛋白质）、库存、电子实验室笔记本和工作流程。

## 何时使用此技能

该技能应该在以下情况下使用：
- 使用 Benchling 的 Python SDK 或 REST API
- 管理生物序列（DNA、RNA、蛋白质）和注册实体
- 自动化库存操作（样品、容器、位置、转移）
- 创建或查询电子实验室笔记本条目
- 构建工作流程自动化或基准测试应用程序
- 在基准测试和外部系统之间同步数据
- 查询基准数据仓库进行分析
- 设置与 AWS EventBridge 的事件驱动集成

## 核心能力

### 1. 身份验证和设置

**Python SDK安装：**
```python
# Stable release
uv pip install benchling-sdk
# or with Poetry
poetry add benchling-sdk
```

**身份验证方法：**

API 密钥身份验证（推荐用于脚本）：
<<<代码块_1>>>

OAuth 客户端凭据（适用于应用程序）：
<<<代码块_2>>>

**要点：**
- API 密钥从 Benchling 中的配置文件设置中获取
- 安全地存储凭据（使用环境变量或密码管理器）
- 所有 API 请求都需要 HTTPS
- 身份验证权限镜像用户界面中的用户权限

有关包括 OIDC 和安全最佳实践在内的详细身份验证信息，请参阅`references/authentication.md`。

### 2. 注册和实体管理

注册实体包括 DNA 序列、RNA 序列、AA 序列、自定义实体和混合物。 SDK 提供了用于创建和管理这些实体的类型化类。

**创建 DNA 序列：**
<<<代码块_3>>>

**注册处注册：**

要在创建时直接注册实体：
<<<代码块_4>>>

**重要提示：** 使用 `entity_registry_id` 或 `naming_strategy`，切勿同时使用。

**更新实体：**
<<<代码块_5>>>

未指定的字段保持不变，允许部分更新。

**列表和分页：**
<<<代码块_6>>>

**关键操作：**
- 创建：`benchling.<entity_type>.create()`
- 读取：`benchling.<entity_type>.get(id)` 或 `.list()`
- 更新：`benchling.<entity_type>.update(id, update_object)`
- 存档：`benchling.<entity_type>.archive(id)`

实体类型：`dna_sequences`、`rna_sequences`、`aa_sequences`、`custom_entities`、`mixtures`

有关全面的 SDK 参考和高级模式，请参阅 `references/sdk_reference.md`。

### 3. 库存管理

管理 Benchling 库存系统中的实物样品、容器、盒子和位置。

**创建容器：**
```python
from benchling_sdk.models import ContainerCreate

container = benchling.containers.create(
    ContainerCreate(
        name="Sample Tube 001",
        schema_id="cont_schema_abc123",
        parent_storage_id="box_abc123",  # optional
        fields=benchling.models.fields({"concentration": "100 ng/μL"})
    )
)
```

**管理盒子：**
```python
from benchling_sdk.models import BoxCreate

box = benchling.boxes.create(
    BoxCreate(
        name="Freezer Box A1",
        schema_id="box_schema_abc123",
        parent_storage_id="loc_abc123"
    )
)
```

**转移物品：**
```python
# Transfer a container to a new location
transfer = benchling.containers.transfer(
    container_id="cont_abc123",
    destination_id="box_xyz789"
)
```

**关键库存操作：**
- 创建容器、盒子、位置、盘子
- 更新库存物品属性
- 在地点之间转移物品
- 签入/签出物品
- 批量传输的批量操作

### 4. 笔记本和文档

与电子实验室笔记本 (ELN) 条目、协议和模板进行交互。

**创建笔记本条目：**
```python
from benchling_sdk.models import EntryCreate

entry = benchling.entries.create(
    EntryCreate(
        name="Experiment 2025-10-20",
        folder_id="fld_abc123",
        schema_id="entry_schema_abc123",
        fields=benchling.models.fields({"objective": "Test gene expression"})
    )
)
```

**将实体链接到条目：**
```python
# Add references to entities in an entry
entry_link = benchling.entry_links.create(
    entry_id="entry_abc123",
    entity_id="seq_xyz789"
)
```

**关键笔记本操作：**
- 创建和更新实验室笔记本条目
- 管理条目模板
- 将实体和结果链接到条目
- 导出文档条目

### 5. 工作流程和自动化

使用 Benchling 的工作流程系统实现实验室流程自动化。

**创建工作流任务：**
```python
from benchling_sdk.models import WorkflowTaskCreate

task = benchling.workflow_tasks.create(
    WorkflowTaskCreate(
        name="PCR Amplification",
        workflow_id="wf_abc123",
        assignee_id="user_abc123",
        fields=benchling.models.fields({"template": "seq_abc123"})
    )
)
```

**更新任务状态：**
```python
from benchling_sdk.models import WorkflowTaskUpdate

updated_task = benchling.workflow_tasks.update(
    task_id="task_abc123",
    workflow_task=WorkflowTaskUpdate(
        status_id="status_complete_abc123"
    )
)
```

**异步操作：**

有些操作是异步的并返回任务：
```python
# Wait for task completion
from benchling_sdk.helpers.tasks import wait_for_task

result = wait_for_task(
    benchling,
    task_id="task_abc123",
    interval_wait_seconds=2,
    max_wait_seconds=300
)
```

**关键工作流程操作：**
- 创建和管理工作流程任务
- 更新任务状态和分配
- 异步执行批量操作
- 监控任务进度

### 6. 事件与集成

使用 AWS EventBridge 订阅基准测试事件以进行实时集成。

**事件类型：**
- 实体创建、更新、归档
- 库存转移
- 工作流任务状态变化
- 条目创建和更新
- 成绩登记

**整合模式：**
1. 在基准设置中配置到 AWS EventBridge 的事件路由
2. 创建EventBridge规则来过滤事件
3. 将事件路由到 Lambda 函数或其他目标
4. 处理事件并更新外部系统

**使用案例：**
- 将基准测试数据同步到外部数据库
- 工作流程完成时触发下游流程
- 发送有关实体更改的通知
- 审计跟踪记录

有关事件架构和配置，请参阅 Benchling 的事件文档。

### 7. 数据仓库和分析

通过数据仓库使用 SQL 查询历史 Benchling 数据。

**访问方式：**
Benchling 数据仓库提供对 Benchling 数据的 SQL 访问以进行分析和报告。使用标准 SQL 客户端和提供的凭据进行连接。

**常见查询：**
- 汇总实验结果
- 分析库存趋势
- 生成合规报告
- 导出数据以供外部分析

**与分析工具集成：**
- 用于交互式分析的 Jupyter 笔记本
- BI 工具（Tableau、Looker、PowerBI）
- 自定义仪表板

## 最佳实践

### 错误处理

SDK会自动重试失败的请求：
```python
# Automatic retry for 429, 502, 503, 504 status codes
# Up to 5 retries with exponential backoff
# Customize retry behavior if needed
from benchling_sdk.retry import RetryStrategy

benchling = Benchling(
    url="https://your-tenant.benchling.com",
    auth_method=ApiKeyAuth("your_api_key"),
    retry_strategy=RetryStrategy(max_retries=3)
)
```

### 分页效率

使用生成器进行内存高效的分页：
```python
# Generator-based iteration
for page in benchling.dna_sequences.list():
    for sequence in page:
        process(sequence)

# Check estimated count without loading all pages
total = benchling.dna_sequences.list().estimated_count()
```

### 模式字段助手

对自定义架构字段使用 `fields()` 帮助器：
```python
# Convert dict to Fields object
custom_fields = benchling.models.fields({
    "concentration": "100 ng/μL",
    "date_prepared": "2025-10-20",
    "notes": "High quality prep"
})
```

### 前向兼容性

SDK 可以优雅地处理未知的枚举值和类型：
- 保留未知的枚举值
- 无法识别的多态类型返回`UnknownType`
- 允许使用较新的 API 版本

### 安全考虑

- 切勿将 API 密钥提交给版本控制
- 使用环境变量作为凭据
- 如果密钥被泄露，请轮换密钥
- 为应用程序授予最低限度的必要权限
- 将 OAuth 用于多用户场景

## 资源

###参考资料/

深入信息的详细参考文档：

- **authentication.md** - 全面的身份验证指南，包括 OIDC、安全最佳实践和凭证管理
- **sdk_reference.md** - 详细的 Python SDK 参考，包含高级模式、示例和所有实体类型
- **api_endpoints.md** - 用于无需 SDK 的直接 HTTP 调用的 REST API 端点参考

根据特定集成要求的需要加载这些参考。

### 脚本/

此技能当前包括示例脚本，可以将其删除或替换为针对特定基准测试工作流程的自定义自动化脚本。

## 常见用例

**1.批量实体导入：**
```python
# Import multiple sequences from FASTA file
from Bio import SeqIO

for record in SeqIO.parse("sequences.fasta", "fasta"):
    benchling.dna_sequences.create(
        DnaSequenceCreate(
            name=record.id,
            bases=str(record.seq),
            is_circular=False,
            folder_id="fld_abc123"
        )
    )
```

**2.库存审核：**
```python
# List all containers in a specific location
containers = benchling.containers.list(
    parent_storage_id="box_abc123"
)

for page in containers:
    for container in page:
        print(f"{container.name}: {container.barcode}")
```

**3.工作流程自动化：**
```python
# Update all pending tasks for a workflow
tasks = benchling.workflow_tasks.list(
    workflow_id="wf_abc123",
    status="pending"
)

for page in tasks:
    for task in page:
        # Perform automated checks
        if auto_validate(task):
            benchling.workflow_tasks.update(
                task_id=task.id,
                workflow_task=WorkflowTaskUpdate(
                    status_id="status_complete"
                )
            )
```

**4.数据导出：**
```python
# Export all sequences with specific properties
sequences = benchling.dna_sequences.list()
export_data = []

for page in sequences:
    for seq in page:
        if seq.schema_id == "target_schema_id":
            export_data.append({
                "id": seq.id,
                "name": seq.name,
                "bases": seq.bases,
                "length": len(seq.bases)
            })

# Save to CSV or database
import csv
with open("sequences.csv", "w") as f:
    writer = csv.DictWriter(f, fieldnames=export_data[0].keys())
    writer.writeheader()
    writer.writerows(export_data)
```

## 其他资源

- **官方文档：** https://docs.benchling.com
- **Python SDK 参考：** https://benchling.com/sdk-docs/
- **API参考：** https://benchling.com/api/reference
- **支持：** [电子邮件受保护]