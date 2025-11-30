<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：工具宇宙
描述：在使用生物信息学、化学信息学、基因组学、结构生物学、蛋白质组学和药物发现的科学研究工具和工作流程时使用此技能。通过此技能，可以访问 600 多种科学工具，包括机器学习模型、数据集、API 和分析包。在搜索科学工具、执行计算生物学工作流程、编写多步骤研究流程、访问 OpenTargets/PubChem/UniProt/PDB/ChEMBL 等数据库、执行研究任务工具发现或将科学计算资源集成到 LLM 工作流程中时使用。
---

# 工具宇宙

## 概述

ToolUniverse 是一个统一的生态系统，通过提供对 600 多种科学资源的标准化访问，使 AI 代理能够充当研究科学家。使用此技能在多个研究领域发现、执行和组合科学工具，包括生物信息学、化学信息学、基因组学、结构生物学、蛋白质组学和药物发现。

**关键能力：**
- 访问 600 多个科学工具、模型、数据集和 API
- 使用自然语言、语义搜索或关键字发现工具
- 通过标准化的人工智能工具交互协议执行工具
- 为复杂的研究问题构建多步骤工作流程
- 通过模型上下文协议 (MCP) 与 Claude 桌面/代码集成

## 何时使用此技能

在以下情况下使用此技能：
- 按功能或领域搜索科学工具（例如，“查找蛋白质结构预测工具”）
- 执行计算生物学工作流程（例如疾病靶标识别、药物发现、基因组学分析）
- 访问科学数据库（OpenTargets、PubChem、UniProt、PDB、ChEMBL、KEGG 等）
- 构建多步骤研究流程（例如，目标发现→结构预测→虚拟筛选）
- 处理生物信息学、化学信息学或结构生物学任务
- 分析基因表达、蛋白质序列、分子结构或临床数据
- 进行文献检索、通路富集或变异注释
- 建立自动化的科学研究工作流程

## 快速入门

### 基本设置
```python
from tooluniverse import ToolUniverse

# Initialize and load tools
tu = ToolUniverse()
tu.load_tools()  # Loads 600+ scientific tools

# Discover tools
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {
        "description": "disease target associations",
        "limit": 10
    }
})

# Execute a tool
result = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {"efoId": "EFO_0000537"}  # Hypertension
})
```

### 模型上下文协议 (MCP)
对于 Claude 桌面/代码集成：
<<<代码块_1>>>

## 核心工作流程

### 1. 工具发现

找到适合您的研究任务的相关工具：

**三种发现方法：**
- `Tool_Finder` - 基于嵌入的语义搜索（需要 GPU）
- `Tool_Finder_LLM` - 基于 LLM 的语义搜索（无需 GPU）
- `Tool_Finder_Keyword` - 快速关键字搜索

**示例：**
<<<代码块_2>>>

**参见`references/tool-discovery.md`了解：**
- 详细的发现方法和搜索策略
- 特定领域的关键词建议
- 寻找工具的最佳实践

### 2. 工具执行

通过标准化接口执行各个工具：

**示例：**
<<<代码块_3>>>

**参见`references/tool-execution.md`了解：**
- 跨领域的真实执行示例
- 工具参数处理和验证
- 结果处理和错误处理
- 生产使用的最佳实践

### 3. 工具组成和工作流程

为复杂的研究工作流程构建多个工具：

**药物发现示例：**
<<<代码块_4>>>

**参见`references/tool-composition.md`了解：**
- 完整的工作流程示例（药物发现、基因组学、临床）
- 顺序和并行工具组合模式
- 输出处理钩子
- 工作流程最佳实践

## 科学领域

ToolUniverse 支持主要科学领域的 600 多种工具：

**生物信息学：**
- 序列分析、比对、BLAST
- 基因表达（RNA-seq、DESeq2）
- 通路富集（KEGG、Reactome、GO）
- 变异注释（VEP、ClinVar）

**化学信息学：**
- 分子描述符和指纹
- 药物发现和虚拟筛选
- ADMET预测和药物相似性
- 化学数据库（PubChem、ChEMBL、ZINC）

**结构生物学：**
- 蛋白质结构预测（AlphaFold）
- 结构检索（PDB）
- 结合位点检测
- 蛋白质-蛋白质相互作用

**蛋白质组学：**
- 质谱分析
- 蛋白质数据库（UniProt、STRING）
- 翻译后修饰

**基因组学：**
- 基因组组装和注释
- 拷贝数变异
- 临床基因组学工作流程

**医学/临床：**
- 疾病数据库（OpenTargets、OMIM）
- 临床试验和 FDA 数据
- 变体分类

**参见`references/domains.md`了解：**
- 完整的域名分类
- 按学科划分的工具示例
- 跨域应用
- 按域搜索策略

## 参考文档

该技能包括全面的参考文件，提供特定方面的详细信息：

- **`references/installation.md`** - 安装、设置、MCP 配置、平台集成
- **`references/tool-discovery.md`>** - 发现方法、搜索策略、列表工具
- **`references/tool-execution.md`** - 执行模式、实际示例、错误处理
- **`references/tool-composition.md`** - 工作流程组合、复杂管道、并行执行
- **`references/domains.md`** - 按领域、用例示例进行工具分类
- **`references/api_reference.md`** - Python API 文档、挂钩、协议

**工作流程：** 在帮助完成特定任务时，请参考相应的文件以获取详细说明。例如，如果搜索工具，请查阅`references/tool-discovery.md`以获取搜索策略。

## 示例脚本

两个可执行示例脚本演示了常见用例：

**`scripts/example_tool_search.py`** - 演示所有三种发现方法：
- 基于关键字的搜索
- 基于法学硕士的搜索
- 特定领域的搜索
- 获取详细的工具信息

**`scripts/example_workflow.py`** - 完整的工作流程示例：
- 药物发现管道（疾病→靶标→结构→筛选→候选药物）
- 基因组分析（表达数据→差异分析→通路）

运行示例以了解典型的使用模式和工作流程组成。

## 最佳实践

1. **工具发现：**
   - 从广泛搜索开始，然后根据结果进行细化
   - 使用 `Tool_Finder_Keyword` 快速搜索已知术语
   - 使用`Tool_Finder_LLM`进行复杂的语义查询
   - 设置适当的`limit`参数（默认值：10）

2. **工具执行：**
   - 执行前始终验证工具参数
   - 为生产工作流程实施错误处理
   - 验证输入数据格式（SMILES、UniProt ID、基因符号）
   - 检查结果类型和结构

3. **工作流程组成：**
   - 在构建完整的工作流程之前单独测试每个步骤
   - 为长工作流程实施检查点
   - 考虑远程 API 的速率限制
   - 当工具独立时使用并行执行

4. **整合：**
   - 初始化 ToolUniverse 一次并重用实例
   - 启动时调用`load_tools()`一次
   - 缓存常用工具信息
   - 启用日志记录以进行调试

## 关键术语

- **工具**：可通过 ToolUniverse 访问的科学资源（模型、数据集、API、包）
- **工具发现**：使用搜索方法（Finder、LLM、关键字）查找相关工具
- **工具执行**：通过 `tu.run()` 运行具有特定参数的工具
- **工具组合**：链接多个工具以实现多步骤工作流程
- **MCP**：模型上下文协议，用于与 Claude 桌面/代码集成
- **AI-Tool交互协议**：LLM-工具通信的标准化接口

## 资源

- **官方网站**：https://aiscientist.tools
- **GitHub**：https://github.com/mims-harvard/ToolUniverse
- **文档**：https://zitniklab.hms.harvard.edu/ToolUniverse/
- **安装**：`uv uv pip install tooluniverse`
- **MCP 服务器**：`tooluniverse-smcp`