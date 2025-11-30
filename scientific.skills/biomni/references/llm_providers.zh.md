<!-- 此文件由机器翻译自 llm_providers.md -->

# LLM 提供商配置

使用 biomni 配置不同 LLM 提供商的综合指南。

## 概述

Biomni 支持多个 LLM 提供商，以便跨不同的基础设施和成本要求进行灵活部署。该框架通过统一的接口抽象出提供商的差异。

## 支持的提供商

1. **人性克劳德**（推荐）
2. **开放人工智能**
3. **Azure OpenAI**
4. **谷歌双子座**
5. **格罗克**
6. **AWS 基岩**
7. **自定义端点**

## 人类克劳德

**推荐用于：** 推理质量、速度和生物医学知识的最佳平衡。

### 设置

```bash
# Set API key
export ANTHROPIC_API_KEY="sk-ant-..."

# Or in .env file
echo "ANTHROPIC_API_KEY=sk-ant-..." >> .env
```

### 可用型号

<<<代码块_1>>>

### 配置选项

<<<代码块_2>>>

**型号特点：**

|型号|最适合 |速度|成本|推理质量|
|--------|----------|--------|------|--------------------|
|作品 4 |复杂的多步骤分析 |慢一点 |高|最高|
|十四行诗 4 |一般生物医学任务|快|中等|高|
|俳句 4 |简单查询，批量处理 |最快|低|好 |

## 开放人工智能

**推荐用于：** 已建立的基础设施、GPT-4 优化。

### 设置

<<<代码块_3>>>

### 可用型号

<<<代码块_4>>>

### 配置

<<<代码块_5>>>

**注意事项：**
- 推荐使用 GPT-4 Turbo 以实现成本效益
- 可能需要额外的生物医学背景来完成专门任务
- 费率限制因账户等级而异

## Azure OpenAI

**推荐用于：** 企业部署、数据驻留要求。

### 设置

<<<代码块_6>>>

### 配置

```python
from biomni.config import default_config

default_config.llm = "azure-gpt-4"
default_config.azure_openai_api_key = "..."
default_config.azure_openai_endpoint = "https://your-resource.openai.azure.com/"
default_config.azure_openai_deployment_name = "gpt-4"
default_config.azure_openai_api_version = "2024-02-01"
```

### 用法

```python
agent = A1(path='./data', llm='azure-gpt-4')
```

**部署注意事项：**
- 需要 Azure OpenAI 服务配置
- 在 Azure 资源创建期间设置的部署名称
- Microsoft定期更新API版本

## 谷歌双子座

**推荐用于：** Google Cloud 集成、多模式任务。

### 设置

```bash
export GOOGLE_API_KEY="..."
```

### 可用型号

```python
# Gemini 2.0 Flash (recommended)
agent = A1(path='./data', llm='gemini-2.0-flash-exp')

# Gemini Pro
agent = A1(path='./data', llm='gemini-pro')
```

### 配置

```python
from biomni.config import default_config

default_config.llm = "gemini-2.0-flash-exp"
default_config.google_api_key = "..."
default_config.llm_temperature = 0.7
```

**特点：**
- 原生多模式支持（文本、图像、代码）
- 快速推理
- 有竞争力的价格

## 格罗克

**推荐用于：** 超快速推理、成本敏感型应用。

### 设置

```bash
export GROQ_API_KEY="gsk_..."
```

### 可用型号

```python
# Llama 3.3 70B
agent = A1(path='./data', llm='llama-3.3-70b-versatile')

# Mixtral 8x7B
agent = A1(path='./data', llm='mixtral-8x7b-32768')
```

### 配置

```python
from biomni.config import default_config

default_config.llm = "llama-3.3-70b-versatile"
default_config.groq_api_key = "gsk_..."
```

**特点：**
- 通过定制硬件进行极快的推理
- 开源模型选项
- 某些模型的上下文窗口有限

## AWS 基岩

**推荐用于：** AWS 基础设施、合规性要求。

### 设置

```bash
export AWS_ACCESS_KEY_ID="..."
export AWS_SECRET_ACCESS_KEY="..."
export AWS_DEFAULT_REGION="us-east-1"
```

### 可用型号

```python
# Claude via Bedrock
agent = A1(path='./data', llm='bedrock-claude-sonnet-4')

# Llama via Bedrock
agent = A1(path='./data', llm='bedrock-llama-3-70b')
```

### 配置

```python
from biomni.config import default_config

default_config.llm = "bedrock-claude-sonnet-4"
default_config.aws_access_key_id = "..."
default_config.aws_secret_access_key = "..."
default_config.aws_region = "us-east-1"
```

**要求：**
- 启用了 Bedrock 访问的 AWS 账户
- 通过AWS控制台请求模型访问
- 为 Bedrock API 配置的 IAM 权限

## 自定义端点

**推荐用于：** 自托管模型、自定义基础设施。

### 配置

```python
from biomni.config import default_config

default_config.llm = "custom"
default_config.custom_llm_endpoint = "http://localhost:8000/v1/chat/completions"
default_config.custom_llm_api_key = "..."  # If required
default_config.custom_llm_model_name = "llama-3-70b"
```

### 用法

```python
agent = A1(path='./data', llm='custom')
```

**端点要求：**
- 必须实现与 OpenAI 兼容的聊天完成 API
- 建议支持函数/工具调用
- JSON 响应格式

**vLLM 示例：**

```bash
# Start vLLM server
python -m vllm.entrypoints.openai.api_server \
    --model meta-llama/Llama-3-70b-chat \
    --port 8000

# Configure biomni
export CUSTOM_LLM_ENDPOINT="http://localhost:8000/v1/chat/completions"
```

## 型号选择指南

### 按任务复杂性

**简单查询**（基因查找、基本计算）：
- 克劳德俳句 4
- 双子座2.0闪存
- 格罗克羊驼 3.3 70B

**中等任务**（数据分析、文献检索）：
- 克劳德十四行诗 4（推荐）
- GPT-4 涡轮增压
- 双子座2.0闪存

**复杂分析**（多步推理，新颖见解）：
- 克劳德作品 4（推荐）
- GPT-4
——克劳德十四行诗 4

### 按成本敏感度

**注重预算：**
1. Groq（最快、最便宜）
2. 克劳德俳句 4
3.双子座2.0闪存

**平衡：**
1.克劳德十四行诗4（推荐）
2.GPT-4涡轮
3.双子座专业版

**质量第一：**
1.克劳德作品4
2.GPT-4
3. 克劳德十四行诗 4

### 按基础设施

**与云无关：**
- 人类克劳德（直接 API）
- OpenAI（直接API）

**AWS 生态系统：**
- AWS 基岩（Claude、Llama）

**Azure 生态系统：**
- Azure OpenAI 服务

**谷歌云：**
- 谷歌双子座

**本地：**
- 具有自托管模型的自定义端点

## 性能比较

基于 Biomni-Eval1 基准：

|供应商|型号|平均分数 |平均时间（秒）|成本/1K 任务 |
|----------|---------|------------|--------------|---------------|
|人择 |作品 4 | 0.89 | 0.89 45 | 45 120 美元 |
|人择 |十四行诗 4 | 0.85 | 0.85 28 | 28 45 美元 |
|开放人工智能 | GPT-4 涡轮 | 0.82 | 0.82 35 | 35 55 美元 |
|谷歌 |双子座2.0闪存| 0.78 | 0.78 22 | 22 25 美元 |
|格罗克 |骆驼 3.3 70B | 0.73 | 0.73 12 | 12 8 美元 |
|人择 |俳句 4 | 0.75 | 0.75 15 | 15 15 美元 |

*注意：成本为近似值，并根据使用模式而变化。*

## 故障排除

### API 密钥问题

```python
# Verify key is set
import os
print(os.getenv('ANTHROPIC_API_KEY'))

# Or check in Python
from biomni.config import default_config
print(default_config.anthropic_api_key)
```

### 速率限制

```python
from biomni.config import default_config

# Add retry logic
default_config.max_retries = 5
default_config.retry_delay = 10  # seconds

# Reduce concurrency
default_config.max_concurrent_requests = 1
```

### 超时错误

```python
# Increase timeout for slow providers
default_config.llm_timeout = 120  # seconds

# Or switch to faster model
default_config.llm = "claude-sonnet-4-20250514"  # Fast and capable
```

### 型号不可用

```bash
# For Bedrock: Enable model access in AWS console
aws bedrock list-foundation-models --region us-east-1

# For Azure: Check deployment name
az cognitiveservices account deployment list \
    --name your-resource-name \
    --resource-group your-rg
```

## 最佳实践

### 成本优化

1. **使用适当的模型** - 不要使用 Opus 4 进行简单查询
2. **启用缓存** - 跨任务重用数据湖访问
3. **批处理** - 将相似的任务分组在一起
4. **监控使用情况** - 跟踪每个任务类型的 API 成本

```python
from biomni.config import default_config

# Enable response caching
default_config.enable_caching = True
default_config.cache_ttl = 3600  # 1 hour
```

### 多提供商策略

```python
def get_agent_for_task(task_complexity):
    """Select provider based on task requirements"""
    if task_complexity == 'simple':
        return A1(path='./data', llm='claude-haiku-4-20250514')
    elif task_complexity == 'moderate':
        return A1(path='./data', llm='claude-sonnet-4-20250514')
    else:
        return A1(path='./data', llm='claude-opus-4-20250514')

# Use appropriate model
agent = get_agent_for_task('moderate')
result = agent.go(task_query)
```

### 后备配置

```python
from biomni.exceptions import LLMError

def execute_with_fallback(task_query):
    """Try multiple providers if primary fails"""
    providers = [
        'claude-sonnet-4-20250514',
        'gpt-4-turbo',
        'gemini-2.0-flash-exp'
    ]

    for llm in providers:
        try:
            agent = A1(path='./data', llm=llm)
            return agent.go(task_query)
        except LLMError as e:
            print(f"{llm} failed: {e}")
            continue

    raise Exception("All providers failed")
```

## 针对提供商的提示

### 人类克劳德
- 最适合复杂的生物医学推理
- 使用 Sonnet 4 完成大多数任务
- 为新颖的研究问题保留 Opus 4

### 开放人工智能
- 添加具有生物医学背景的系统提示以获得更好的结果
- 使用 JSON 模式进行结构化输出
- 监控令牌使用情况 - 上下文窗口限制

### Azure OpenAI
- 在靠近数据的区域提供部署
- 使用托管身份进行安全身份验证
- 监视 Azure 门户中的配额消耗

### 谷歌双子座
- 利用多模式功能执行基于图像的任务
- 使用流式传输进行长时间运行的分析
- 考虑使用 Gemini Pro 来处理生产工作负载

### 格罗克
- 高通量筛选任务的理想选择
- 与 Claude/GPT-4 相比推理深度有限
- 最适合定义明确的结构化问题

### AWS 基岩
- 尽可能使用 IAM 角色而不是访问密钥
- 启用 CloudWatch 日志记录以进行调试
- 监控跨区域延迟