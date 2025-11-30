<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：困惑搜索
描述：通过 LiteLLM 和 OpenRouter 使用 Perplexity 模型执行基于 AI 的网络搜索和实时信息。在对当前信息进行网络搜索、查找最新的科学文献、通过来源引用获得有根据的答案或访问超出模型知识范围的信息时，应该使用此技能。通过单个 OpenRouter API 密钥提供对多个 Perplexity 模型的访问，包括 Sonar Pro、Sonar Pro Search（高级代理搜索）和 Sonar Reasoning Pro。
---

# 困惑度搜索

## 概述

通过 LiteLLM 和 OpenRouter 使用 Perplexity 模型执行人工智能驱动的网络搜索。 Perplexity 提供实时、基于网络的答案以及来源引文，使其成为查找当前信息、最新科学文献以及超出模型训练数据截止范围的事实的理想选择。

此技能可以通过 OpenRouter 访问所有 Perplexity 模型，仅需要一个 API 密钥（不需要单独的 Perplexity 帐户）。

## 何时使用此技能

在以下情况下使用此技能：
- 搜索当前信息或最新动态（2024 年及以后）
- 查找最新的科学出版物和研究
- 获取基于网络资源的实时答案
- 通过来源引用验证事实
- 跨多个领域进行文献检索
- 访问超出模型知识范围的信息
- 进行特定领域的研究（生物医学、技术、临床）
- 比较当前的方法或技术

**请勿将**用于：
- 简单的计算或逻辑问题（直接使用）
- 需要代码执行的任务（使用标准工具）
- 模型训练数据范围内的问题（除非需要验证）

## 快速入门

### 设置（一次性）

1. **获取 OpenRouter API 密钥**：
   - 访问https://openrouter.ai/keys
   - 创建帐户并生成API密钥
   - 添加积分到帐户（建议至少 5 美元）

2. **配置环境**：
   ```bash
   # Set API key
   export OPENROUTER_API_KEY='sk-or-v1-your-key-here'

   # Or use setup script
   python scripts/setup_env.py --api-key sk-or-v1-your-key-here
   ```

3. **安装依赖**：
   <<<代码块_1>>>

4. **验证设置**：
   <<<代码块_2>>>

请参阅 `references/openrouter_setup.md` 了解详细的设置说明、故障排除和安全最佳实践。

### 基本用法

**简单搜索：**
<<<代码块_3>>>

**保存结果：**
<<<代码块_4>>>

**使用具体型号：**
<<<代码块_5>>>

**详细输出：**
<<<代码块_6>>>

## 可用型号

通过 `--model` 参数访问模型：

- **sonar-pro**（默认）：通用搜索，成本和质量的最佳平衡
- **sonar-pro-search**：最先进的代理搜索，具有多步骤推理
- **声纳**：基本模型，对于简单查询最具成本效益
- **sonar-reasoning-pro**：通过逐步分析进行高级推理
- **声纳推理**：基本推理能力

**选型指南：**
- 默认查询 → `sonar-pro`
- 复杂的多步骤分析 → `sonar-pro-search`
- 需要明确的推理 → `sonar-reasoning-pro`
- 简单的事实查找 → `sonar`
- 成本敏感的批量查询 → `sonar`

有关详细比较、用例、定价和性能特征，请参阅`references/model_comparison.md`。

## 制定有效的查询

### 具体而详细

**好的例子：**
- “2024年发表的CAR-T细胞疗法治疗B细胞淋巴瘤的最新临床试验结果是什么？”
-“比较 mRNA 疫苗与病毒载体疫苗针对 COVID-19 的功效和安全性”
-“通过 2023-2024 年研究的具体准确度指标解释 AlphaFold3 相对于 AlphaFold2 的改进”

**不好的例子：**
- “告诉我有关癌症治疗的信息”（太宽泛）
-“CRISPR”（太模糊）
- “疫苗”（缺乏特异性）

### 包括时间限制

Perplexity 搜索实时网络数据：
- “2024 年《自然医学》发表了哪些关于长期新冠肺炎的论文？”
- “大型语言模型效率的最新进展（过去 6 个月）是什么？”
- “NeurIPS 2023 上宣布了哪些关于人工智能安全的内容？”

### 指定域和来源

为了获得高质量的结果，请提及来源偏好：
- “根据高影响力期刊上经过同行评审的出版物......”
- “基于 FDA 批准的治疗方法......”
- “来自临床试验注册中心，例如 ClinicalTrials.gov...”

### 构造复杂查询

将复杂的问题分解为清晰的部分：
1. **主题**：主要主题
2. **范围**：感兴趣的具体方面
3. **背景**：时间框架、领域、限制
4. **输出**：所需的答案格式或类型
**示例：**
“根据 2023 年至 2024 年发表的研究，AlphaFold3 相对于 AlphaFold2 在蛋白质结构预测方面有何改进？包括具体的准确性指标和基准。”

有关查询设计、特定领域模式和高级技术的综合指南，请参阅 `references/search_strategies.md`。

## 常见用例

### 科学文献检索

```bash
python scripts/perplexity_search.py \
  "What does recent research (2023-2024) say about the role of gut microbiome in Parkinson's disease? Focus on peer-reviewed studies and include specific bacterial species identified." \
  --model sonar-pro
```

### 技术文档

```bash
python scripts/perplexity_search.py \
  "How to implement real-time data streaming from Kafka to PostgreSQL using Python? Include considerations for handling backpressure and ensuring exactly-once semantics." \
  --model sonar-reasoning-pro
```

### 对比分析

```bash
python scripts/perplexity_search.py \
  "Compare PyTorch versus TensorFlow for implementing transformer models in terms of ease of use, performance, and ecosystem support. Include benchmarks from recent studies." \
  --model sonar-pro-search
```

### 临床研究

```bash
python scripts/perplexity_search.py \
  "What is the evidence for intermittent fasting in managing type 2 diabetes in adults? Focus on randomized controlled trials and report HbA1c changes and weight loss outcomes." \
  --model sonar-pro
```

### 趋势分析

```bash
python scripts/perplexity_search.py \
  "What are the key trends in single-cell RNA sequencing technology over the past 5 years? Highlight improvements in throughput, cost, and resolution, with specific examples." \
  --model sonar-pro
```

## 处理结果

### 编程访问

使用 `perplexity_search.py` 作为模块：

```python
from scripts.perplexity_search import search_with_perplexity

result = search_with_perplexity(
    query="What are the latest CRISPR developments?",
    model="openrouter/perplexity/sonar-pro",
    max_tokens=4000,
    temperature=0.2,
    verbose=False
)

if result["success"]:
    print(result["answer"])
    print(f"Tokens used: {result['usage']['total_tokens']}")
else:
    print(f"Error: {result['error']}")
```

### 保存并处理结果

```bash
# Save to JSON
python scripts/perplexity_search.py "query" --output results.json

# Process with jq
cat results.json | jq '.answer'
cat results.json | jq '.usage'
```

### 批处理

创建用于多个查询的脚本：

```bash
#!/bin/bash
queries=(
  "CRISPR developments 2024"
  "mRNA vaccine technology advances"
  "AlphaFold3 accuracy improvements"
)

for query in "${queries[@]}"; do
  echo "Searching: $query"
  python scripts/perplexity_search.py "$query" --output "results_$(echo $query | tr ' ' '_').json"
  sleep 2  # Rate limiting
done
```

## 成本管理

困惑模型有不同的定价等级：

**每次查询的大约成本：**
- 声纳：0.001-0.002美元（最具成本效益）
- Sonar Pro：0.002-0.005 美元（建议默认）
- 声纳推理专业版：0.005-0.010 美元
- Sonar Pro 搜索：$0.020-0.050+（最全面）

**成本优化策略：**
1. 使用 `sonar` 进行简单的事实查找
2. 对于大多数查询，默认为 `sonar-pro`
3. 保留`sonar-pro-search`用于复杂分析
4. 设置`--max-tokens`限制响应长度
5. 监控 https://openrouter.ai/activity 的使用情况
6. 在 OpenRouter 仪表板中设置支出限额

## 故障排除

### API 密钥未设置

**错误**：“未配置 OpenRouter API 密钥”

**解决方案**：
```bash
export OPENROUTER_API_KEY='sk-or-v1-your-key-here'
# Or run setup script
python scripts/setup_env.py --api-key sk-or-v1-your-key-here
```

### LiteLLM 未安装

**错误**：“未安装 LiteLLM”

**解决方案**：
```bash
uv pip install litellm
```

### 速率限制

**错误**：“超出速率限制”

**解决方案**：
- 等待几秒钟再重试
- 提高 https://openrouter.ai/keys 的速率限制
- 在批处理中的请求之间添加延迟

### 学分不足

**错误**：“积分不足”

**解决方案**：
- 在 https://openrouter.ai/account 添加积分
- 启用自动充电以防止中断

请参阅 `references/openrouter_setup.md` 获取全面的故障排除指南。

## 与其他技能的整合

该技能补充了其他科学技能：

###文献综述

与 `literature-review` 技能一起使用：
1. 使用 Perplexity 查找最近的论文和预印本
2. 用实时网络结果补充 PubMed 搜索
3. 验证引用并查找相关工作
4. 发现数据库索引后的最新进展

### 科学写作

与 `scientific-writing` 技能一起使用：
1. 查找最近的参考文献以进行介绍/讨论
2. 验证当前的技术水平
3.检查最新术语和约定
4. 确定最近的竞争方法

### 假设生成

与 `hypothesis-generation` 技能一起使用：
1. 搜索最新研究成果
2. 确定当前的知识差距
3. 发现最新的方法论进展
4.发现新兴研究方向

### 批判性思维

与 `scientific-critical-thinking` 技能一起使用：
1. 寻找支持和反对假设的证据
2. 找到方法论批评
3. 识别领域内的争议
4. 用现有证据验证主张

## 最佳实践

### 查询设计

1. **具体**：包括领域、时间范围和限制
2. **使用术语**：适合领域的关键词和短语
3. **指定来源**：提及首选的出版物类型或期刊
4. **结构问题**：具有明确上下文的清晰组件
5. **迭代**：根据初始结果进行细化

### 型号选择

1. **从 sonar-pro 开始**：对于大多数查询来说，默认值很好
2. **复杂性升级**：使用sonar-pro-search进行多步分析
3. **为简单起见降级**：使用声纳了解基本事实
4. **使用推理模型**：当需要逐步分析时

### 成本优化

1. **选择合适的模型**：将模型与查询复杂度相匹配
2. **设置代币限制**：使用`--max-tokens`来控制成本
3. **监控使用情况**：定期检查OpenRouter仪表板
4. **高效批处理**：尽可能组合相关的简单查询
5. **缓存结果**：保存并重用重复查询的结果

### 安全

1. **保护 API 密钥**：永远不要进行版本控制
2. **使用环境变量**：将密钥与代码分开
3. **设置支出限额**：在 OpenRouter 仪表板中配置
4. **监控使用情况**：注意意外活动
5. **轮换密钥**：定期更换密钥

## 资源

### 捆绑资源

**脚本：**
- `scripts/perplexity_search.py`：带有 CLI 界面的主搜索脚本
- `scripts/setup_env.py`：环境设置和验证助手

**参考资料：**
- `references/search_strategies.md`：综合查询设计指南
- `references/model_comparison.md`：详细的型号比较和选型指南
- `references/openrouter_setup.md`：完整的设置、故障排除和安全指南

**资产：**
- `assets/.env.example`：示例环境文件模板

### 外部资源

**开放路由器：**
- 仪表板：https://openrouter.ai/account
- API 密钥：https://openrouter.ai/keys
- 困惑模型：https://openrouter.ai/perplexity
- 使用情况监控：https://openrouter.ai/activity
- 文档：https://openrouter.ai/docs

**精简法学硕士：**
- 文档：https://docs.litellm.ai/
- OpenRouter 提供商：https://docs.litellm.ai/docs/providers/openrouter
- GitHub：https://github.com/BerriAI/litellm

**困惑：**
- 官方文档：https://docs.perplexity.ai/

## 依赖关系

### 必填

```bash
# LiteLLM for API access
uv pip install litellm
```

### 可选

```bash
# For .env file support
uv pip install python-dotenv

# For JSON processing (usually pre-installed)
uv pip install jq
```

### 环境变量

要求：
- `OPENROUTER_API_KEY`：您的 OpenRouter API 密钥

可选：
- `DEFAULT_MODEL`：使用的默认模型（默认：sonar-pro）
- `DEFAULT_MAX_TOKENS`：默认最大令牌（默认：4000）
- `DEFAULT_TEMPERATURE`：默认温度（默认值：0.2）

## 总结

该技能提供：

1. **实时网络搜索**：访问训练数据截止之外的当前信息
2. **多种型号**：从经济高效的Sonar到先进的Sonar Pro Search
3. **设置简单**：单个 OpenRouter API 密钥，无需单独的 Perplexity 帐户
4. **全面指导**：查询设计和模型选择的详细参考
5. **成本效益**：即用即付定价，并具有使用情况监控
6. **科学焦点**：针对研究、文献检索和技术查询进行优化
7. **轻松集成**：与其他科学技能无缝协作

进行人工智能驱动的网络搜索，以查找当前信息、最新研究以及带有来源引文的有根据的答案。