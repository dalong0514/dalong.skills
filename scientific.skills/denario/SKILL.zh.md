<!-- 此文件由机器翻译自 SKILL.md -->

---
姓名：迪纳里奥
描述：用于科学研究协助的多智能体人工智能系统，可自动化从数据分析到出版的研究工作流程。当从数据集产生研究想法、开发研究方法、执行计算实验、进行文献检索或以 LaTeX 格式生成可发表的论文时，应该使用此技能。通过可定制的代理编排支持端到端研究管道。
---

# 德纳里奥

## 概述

Denario 是一个多智能体人工智能系统，旨在自动化科学研究工作流程，从初始数据分析到可发表的手稿。它基于 AG2 和 LangGraph 框架构建，协调多个专门代理来处理假设生成、方法开发、计算分析和论文写作。

## 何时使用此技能

在以下情况下使用此技能：
- 分析数据集以生成新颖的研究假设
- 开发结构化研究方法
- 执行计算实验并生成可视化
- 进行文献检索以了解研究背景
- 根据研究成果撰写期刊格式的 LaTeX 论文
- 自动化从数据到发表的完整研究流程

## 安装

使用 uv 安装 denario（推荐）：

```bash
uv init
uv add "denario[app]"
```

或者使用点：

<<<代码块_1>>>

对于 Docker 部署或从源代码构建，请参阅 `references/installation.md`。

## LLM API 配置

Denario 需要来自受支持的 LLM 提供商的 API 密钥。支持的提供商包括：
- 谷歌顶点人工智能
- 开放人工智能
- 与AG2/LangGraph兼容的其他LLM服务

使用环境变量或 `.env` 文件安全地存储 API 密钥。有关包括 Vertex AI 设置在内的详细配置说明，请参阅 `references/llm_configuration.md`。

## 核心研究工作流程

Denario 遵循结构化的四阶段研究流程：

### 1. 数据说明

通过指定可用数据和工具来定义研究背景：

<<<代码块_2>>>

### 2.创意产生

根据数据描述生成研究假设：

<<<代码块_3>>>

这会根据所描述的数据产生研究问题或假设。或者，提供自定义想法：

<<<代码块_4>>>

### 3. 方法开发

制定研究方法：

<<<代码块_5>>>

这创建了一种用于研究假设的结构化方法。还可以接受具有自定义方法的 Markdown 文件：

<<<代码块_6>>>

### 4. 结果生成

执行计算实验并生成分析：

```python
den.get_results()
```

它运行方法、执行计算、创建可视化并生成结果。还可以提供预先计算的结果：

```python
den.set_results("path/to/results.md")
```

### 5. 论文生成

创建可供发表的 LaTeX 论文：

```python
from denario import Journal

den.get_paper(journal=Journal.APS)
```

生成的论文包括指定期刊的正确格式、集成的图形和完整的 LaTeX 源代码。

## 可用期刊

Denario 支持多种日志格式样式：
- `Journal.APS` - 美国物理学会格式
- 可能会提供其他期刊；检查 `references/research_pipeline.md` 以获得完整列表

## 启动图形用户界面

运行图形用户界面：

```bash
denario run
```

这启动了一个基于网络的界面，用于交互式研究工作流程管理。

## 常见工作流程

### 端到端研究管道

```python
from denario import Denario, Journal

# Initialize project
den = Denario(project_dir="./research_project")

# Define research context
den.set_data_description("""
Dataset: Time-series measurements of [phenomenon]
Available tools: pandas, sklearn, scipy
Research goal: Investigate [research question]
""")

# Generate research idea
den.get_idea()

# Develop methodology
den.get_method()

# Execute analysis
den.get_results()

# Create publication
den.get_paper(journal=Journal.APS)
```

### 混合工作流程（自定义+自动化）

```python
# Provide custom research idea
den.set_idea("Investigate the correlation between X and Y using time-series analysis")

# Auto-generate methodology
den.get_method()

# Auto-generate results
den.get_results()

# Generate paper
den.get_paper(journal=Journal.APS)
```

### 文献检索整合

有关文献搜索功能和其他工作流程示例，请参阅`references/examples.md`。

## 高级功能

- **多代理编排**：AG2和LangGraph协调专门的代理来完成不同的研究任务
- **可重复的研究**：所有阶段都会产生可以进行版本控制的结构化输出
- **期刊集成**：目标出版地点的自动格式化
- **灵活输入**：每个管道阶段手动或自动
- **Docker 部署**：具有 LaTeX 和所有依赖项的容器化环境

## 详细参考资料

对于综合文档：
- **安装选项**：`references/installation.md`
- **LLM配置**：`references/llm_configuration.md`
- **完整 API 参考**：`references/research_pipeline.md`
- **工作流程示例**：`references/examples.md`

## 故障排除

常见问题及解决方案：
- **API 密钥错误**：确保环境变量设置正确（参见 `references/llm_configuration.md`）
- **LaTeX 编译**：安装 TeX 发行版或使用预装 LaTeX 的 Docker 镜像
- **包冲突**：使用虚拟环境或Docker进行隔离
- **Python版本**：需要Python 3.12或更高版本