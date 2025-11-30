<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：latchbio整合
描述：“用于生物信息学工作流程的 Latch 平台。使用 Latch SDK、@workflow/@task 装饰器构建管道，部署无服务器工作流程、LatchFile/LatchDir、Nextflow/Snakemake 集成。”
---

# LatchBio 集成

## 概述

Latch 是一个 Python 框架，用于将生物信息学工作流程构建和部署为无服务器管道。基于 Flyte 构建，使用 @workflow/@task 装饰器创建工作流，使用 LatchFile/LatchDir 管理云数据，配置资源并集成 Nextflow/Snakemake 管道。

## 核心能力

Latch 平台提供四个主要功能领域：

### 1. 工作流程创建和部署
- 使用 Python 装饰器定义无服务器工作流程
- 支持原生 Python、Nextflow 和 Snakemake 管道
- 使用 Docker 自动容器化
- 自动生成的无代码用户界面
- 版本控制和可重复性

### 2. 数据管理
- 云存储抽象（LatchFile、LatchDir）
- 具有注册表的结构化数据组织（项目→表格→记录）
- 使用链接和枚举进行类型安全的数据操作
- 本地与云端文件自动传输
- 用于文件选择的全局模式匹配

### 3.资源配置
- 预配置的任务装饰器（@small_task、@large_task、@small_gpu_task、@large_gpu_task）
- 自定义资源规格（CPU、内存、GPU、存储）
- GPU 支持（K80、V100、A100）
- 超时和存储配置
- 成本优化策略

### 4. 验证的工作流程
- 生产就绪的预建管道
- Bulk RNA-seq、DESeq2、通路分析
- 用于蛋白质结构预测的 AlphaFold 和 ColabFold
- 单细胞工具（ArchR、scVelo、emptyDropsR）
- CRISPR 分析、系统发育学等

## 快速入门

### 安装和设置

```bash
# Install Latch SDK
python3 -m uv pip install latch

# Login to Latch
latch login

# Initialize a new workflow
latch init my-workflow

# Register workflow to platform
latch register my-workflow
```

**先决条件：**
- Docker已安装并运行
- 锁存帐户凭据
-Python 3.8+

### 基本工作流程示例

<<<代码块_1>>>

## 何时使用此技能

当遇到以下任一场景时，应使用此技能：

**工作流程开发：**
- “创建用于 RNA-seq 分析的 Latch 工作流程”
- “将我的管道部署到 Latch”
- “将我的 Nextflow 管道转换为 Latch”
- “将 GPU 支持添加到我的工作流程中”
- 使用 `@workflow`、`@task` 装饰器

**数据管理：**
- “在 Latch 注册表中组织我的测序数据”
- “如何使用 LatchFile 和 LatchDir？”
- “在 Latch 中设置样本跟踪”
- 使用 `latch:///` 路径

**资源配置：**
- “为 Latch 上的 AlphaFold 配置 GPU”
- “我的任务内存不足”
- “如何优化工作流程成本？”
- 使用任务装饰器

**经过验证的工作流程：**
- “在闩锁上运行 AlphaFold”
- “使用DESeq2进行差异表达”
- “可用的预建工作流程”
- 使用`latch.verified`模块

## 详细文档

该技能包括按能力组织的综合参考文档：

### 参考文献/workflow-creation.md
**阅读本文的目的是：**
- 创建和注册工作流程
- 任务定义和装饰器
- 支持Python、Nextflow、Snakemake
- 启动计划和有条件的部分
- 工作流程执行（CLI 和编程）
- 多步骤和并行管道
- 解决注册问题

**关键主题：**
- `latch init` 和 `latch register` 命令
- `@workflow` 和 `@task` 装饰器
- LatchFile 和 LatchDir 基础知识
- 类型注释和文档字符串
- 具有预设参数的启动计划
- 有条件的 UI 部分

### 参考文献/data-management.md
**阅读本文的目的是：**
- 带有 LatchFile 和 LatchDir 的云存储
- 登记系统（项目、表格、记录）
- 链接的记录和关系
- 枚举和类型列
- 批量操作和交易
- 与工作流程集成
- 帐户和工作区管理

**关键主题：**
- `latch:///` 路径格式
- 文件传输和全局模式
- 创建和查询注册表
- 列类型（字符串、数字、文件、链接、枚举）
- 记录CRUD操作
- 工作流程-注册表集成

### 参考文献/resource-configuration.md
**阅读本文的目的是：**
- 任务资源装饰器
- 自定义CPU、内存、GPU配置
- GPU 类型（K80、V100、A100）
- 超时和存储设置
- 资源优化策略
- 具有成本效益的工作流程设计
- 监控和调试

**关键主题：**
- `@small_task`、`@large_task`、`@small_gpu_task`、`@large_gpu_task`
- `@custom_task` 具有精确的规格
- 多GPU配置
- 按工作负载类型选择资源
- 平台限制和配额

### 参考文献/verified-workflows.md
**阅读本文的目的是：**
- 预建的生产工作流程
- 批量 RNA 测序和 DESeq2
- AlphaFold 和 ColabFold
- 单细胞分析（ArchR、scVelo）
- CRISPR编辑分析
- 途径丰富
- 与自定义工作流程集成

**关键主题：**
- `latch.verified` 模块导入
- 可用的经过验证的工作流程
- 工作流程参数和选项
- 结合验证和自定义步骤
- 版本管理

## 常见工作流程模式

### 完整的 RNA-seq 流程

<<<代码块_2>>>

### GPU 加速的工作流程

<<<代码块_3>>>

### 注册表集成工作流程

<<<代码块_4>>>

## 最佳实践

### 工作流程设计
1.对所有参数使用类型注释
2. 编写清晰的文档字符串（出现在 UI 中）
3. 从标准任务装饰器开始，根据需要进行扩展
4. 将复杂的工作流程分解为模块化任务
5. 实施适当的错误处理

### 数据管理
6.使用一致的文件夹结构
7. 在批量输入之前定义注册表架构
8. 使用链接记录建立关系
9. 将元数据存储在注册表中以供追溯

### 资源配置
10.适当大小的资源（不要过度分配）
11.仅在算法支持时才使用GPU
12. 监控执行指标并优化
13. 尽可能设计并行执行

### 开发工作流程
14. 注册前使用 Docker 进行本地测试
15.对工作流代码使用版本控制
16. 记录资源需求
17. 分析工作流程以确定实际需求

## 故障排除

### 常见问题

**注册失败：**
- 确保 Docker 正在运行
- 使用`latch login`检查身份验证
- 验证 Dockerfile 中的所有依赖项
- 使用 `--verbose` 标志来获取详细日志

**资源问题：**
- 内存不足：增加任务装饰器中的内存
- 超时：增加超时参数
- 存储问题：增加临时存储_gib

**数据访问：**
- 使用正确的`latch:///`路径格式
- 验证工作区中是否存在文件
- 检查共享工作区的权限

**类型错误：**
- 为所有参数添加类型注释
- 使用 LatchFile/LatchDir 作为文件/目录参数
- 确保工作流程返回类型与实际返回匹配

## 其他资源

- **官方文档**：https://docs.latch.bio
- **GitHub 存储库**：https://github.com/latchbio/latch
- **Slack 社区**：加入 Latch SDK 工作区
- **API 参考**：https://docs.latch.bio/api/latch.html
- **博客**：https://blog.latch.bio

## 支持

对于问题或疑问：
1.检查上面的文档链接
2.搜索GitHub问题
3. 在 Slack 社区提问
4.联系support@latch.bio