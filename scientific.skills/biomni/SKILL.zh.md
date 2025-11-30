<!-- 此文件由机器翻译自 SKILL.md -->

---
姓名：比奥姆尼
描述：自主生物医学人工智能代理框架，用于执行基因组学、药物发现、分子生物学和临床分析等复杂的研究任务。在进行多步骤生物医学研究时使用此技能，包括 CRISPR 筛选设计、单细胞 RNA-seq 分析、ADMET 预测、GWAS 解释、罕见疾病诊断或实验室方案优化。利用法学硕士推理与代码执行和集成生物医学数据库。
---

#比奥姆尼

## 概述

Biomni 是斯坦福大学 SNAP 实验室的开源生物医学 AI 代理框架，可自主执行跨生物医学领域的复杂研究任务。在处理多步骤生物推理任务、分析生物医学数据或进行基因组学、药物发现、分子生物学和临床分析研究时，可以使用此技能。

## 核心能力

Biomni 擅长：

1. **多步生物推理** - 复杂生物医学查询的自主任务分解和规划
2. **代码生成和执行** - 动态分析管道创建以进行数据处理
3. **知识检索** - 访问约 11GB 的综合生物医学数据库和文献
4. **跨领域问题解决** - 基因组学、蛋白质组学、药物发现和临床任务的统一界面

## 何时使用此技能

使用 Biomni 用于：
- **CRISPR筛选** - 设计筛选、优先考虑基因、分析敲除效应
- **单细胞 RNA-seq** - 细胞类型注释、差异表达、轨迹分析
- **药物发现** - ADMET 预测、靶点识别、化合物优化
- **GWAS 分析** - 变异解释、因果基因识别、通路富集
- **临床基因组学** - 罕见疾病诊断、变异致病性、表型-基因型作图
- **实验室协议** - 协议优化、文献综合、实验设计

## 快速入门

### 安装和设置

安装 Biomni 并为 LLM 提供商配置 API 密钥：

```bash
uv pip install biomni --upgrade
```

配置 API 密钥（存储在 `.env` 文件或环境变量中）：
<<<代码块_1>>>

使用 `scripts/setup_environment.py` 获取交互式设置帮助。

### 基本使用模式

<<<代码块_2>>>

## 与 Biomni 合作

### 1.代理初始化

A1 类是 biomni 的主要接口：

<<<代码块_3>>>

**支持的法学硕士提供者：**
- 人类克劳德（推荐）：`claude-sonnet-4-20250514`、`claude-opus-4-20250514`
- OpenAI：`gpt-4`、`gpt-4-turbo`
- Azure OpenAI：通过 Azure 配置
- 谷歌双子座：`gemini-2.0-flash-exp`
- Groq：`llama-3.3-70b-versatile`
- AWS Bedrock：通过 Bedrock API 的各种模型

请参阅 `references/llm_providers.md` 了解详细的 LLM 配置说明。

### 2. 任务执行工作流程

Biomni 遵循自主代理工作流程：

<<<代码块_4>>>

### 3. 常见任务模式

#### CRISPR 筛选设计
<<<代码块_5>>>

#### 单细胞 RNA-seq 分析
<<<代码块_6>>>

#### 药物 ADMET 预测
```python
agent.go("""
Predict ADMET properties for these drug candidates:
[SMILES strings or compound IDs]
Focus on:
- Absorption (Caco-2 permeability, HIA)
- Distribution (plasma protein binding, BBB penetration)
- Metabolism (CYP450 interaction)
- Excretion (clearance)
- Toxicity (hERG liability, hepatotoxicity)
""")
```

#### GWAS 变异解释
```python
agent.go("""
Interpret GWAS results for [trait/disease]:
- Identify genome-wide significant variants
- Map variants to causal genes
- Perform pathway enrichment analysis
- Predict functional consequences
Summary statistics file: [path/to/gwas_summary.txt]
""")
```

有关所有生物医学领域的综合任务示例，请参阅 `references/use_cases.md`。

### 4. 数据集成

Biomni 集成了约 11GB 的生物医学知识源：
- **基因数据库** - Ensembl、NCBI Gene、UniProt
- **蛋白质结构** - PDB、AlphaFold
- **临床数据集** - ClinVar、OMIM、HPO
- **文献索引** - PubMed 摘要、生物医学本体论
- **通路数据库** - KEGG、Reactome、GO

首次使用时，数据会自动下载到指定的 `path`。

### 5.MCP 服务器集成

通过模型上下文协议使用外部工具扩展 biomni：

```python
# MCP servers can provide:
# - FDA drug databases
# - Web search for literature
# - Custom biomedical APIs
# - Laboratory equipment interfaces

# Configure MCP servers in .biomni/mcp_config.json
```

### 6. 评估框架

生物医学任务上代理性能的基准测试：

```python
from biomni.eval import BiomniEval1

evaluator = BiomniEval1()

# Evaluate on specific task types
score = evaluator.evaluate(
    task_type='crispr_design',
    instance_id='test_001',
    answer=agent_output
)

# Access evaluation dataset
dataset = evaluator.load_dataset()
```

## 最佳实践

### 任务制定
- **具体** - 包括生物背景、生物体、细胞类型、条件
- **指定输出** - 明确说明所需的分析输出和格式
- **提供数据路径** - 包括要分析的数据集的文件路径
- **设置限制** - 提及时间/计算限制（如果相关）

### 安全考虑
⚠️ **重要**：Biomni 使用完整的系统权限执行 LLM 生成的代码。对于生产用途：
- 在隔离环境（Docker、VM）中运行
- 避免暴露敏感凭证
- 在敏感上下文中执行之前检查生成的代码
- 尽可能使用沙盒执行环境

### 性能优化
- **选择合适的法学硕士** - 建议使用 Claude Sonnet 4 以平衡速度/质量
- **设置合理的超时** - 调整复杂任务的`default_config.timeout_seconds`
- **监控迭代** - 跟踪 `max_iterations` 以防止失控循环
- **缓存数据** - 跨会话重用下载的数据湖

### 结果文档
```python
# Always save conversation history for reproducibility
agent.save_conversation_history("results/project_name_YYYYMMDD.pdf")

# Include in reports:
# - Original task description
# - Generated analysis code
# - Results and interpretations
# - Data sources used
```

## 资源

### 参考文献
详细文档可在 `references/` 目录中找到：

- **`api_reference.md`** - A1 类、配置和评估的完整 API 文档
- **`llm_providers.md`** - LLM 提供商设置（Anthropic、OpenAI、Azure、Google、Groq、AWS）
- **`use_cases.md`** - 所有生物医学领域的综合任务示例

### 脚本
`scripts/` 目录中的帮助程序脚本：

- **`setup_environment.py`** - 交互环境和 API 密钥配置
- **`generate_report.py`** - 使用自定义格式增强 PDF 报告生成

### 外部资源
- **GitHub**：https://github.com/snap-stanford/biomni
- **网络平台**：https://biomni.stanford.edu
- **论文**：https://www.biorxiv.org/content/10.1101/2025.05.30.656746v1
- **型号**：https://huggingface.co/biomni/Biomni-R0-32B-Preview
- **评估数据集**：https://huggingface.co/datasets/biomni/Eval1

## 故障排除

### 常见问题

**数据下载失败**
```python
# Manually trigger data lake download
agent = A1(path='./data', llm='your-llm')
# First .go() call will download data
```

**API 密钥错误**
```bash
# Verify environment variables
echo $ANTHROPIC_API_KEY
# Or check .env file in working directory
```

**复杂任务超时**
```python
from biomni.config import default_config
default_config.timeout_seconds = 3600  # 1 hour
```

**大型数据集的内存问题**
- 对大文件使用流式传输
- 分块处理数据
- 增加系统内存分配

### 获取帮助

对于问题或疑问：
- GitHub 问题：https://github.com/snap-stanford/biomni/issues
- 文档：检查 `references/` 文件以获取详细指导
- 社区：斯坦福 SNAP 实验室和 biomni 贡献者