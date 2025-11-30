<!-- 此文件由机器翻译自 openrouter_setup.md -->

# OpenRouter 设置指南

设置和使用 OpenRouter 进行 Perplexity 模型访问的完整指南。

## 什么是 OpenRouter？

OpenRouter 是一个统一的 API 网关，可通过单个 API 接口访问来自不同提供商的 100 多个 AI 模型。它提供：

- **单一 API 密钥**：一键访问多个模型
- **统一格式**：OpenAI兼容的API格式
- **成本跟踪**：内置使用情况监控和计费
- **模型路由**：智能回退和负载平衡
- **即用即付**：无需订阅，只需按使用量付费

特别是对于 Perplexity 模型，OpenRouter 提供对某些模型（例如 Sonar Pro Search）的独占访问。

## 开始使用

### 第 1 步：创建 OpenRouter 帐户

1.访问https://openrouter.ai/
2.点击右上角“注册”
3. 使用 Google、GitHub 或电子邮件注册
4. 如果使用电子邮件注册，请验证您的电子邮件

### 第 2 步：添加付款方式

OpenRouter 采用按量付费计费方式：

1. 导航至https://openrouter.ai/account
2. 单击“制作人员”选项卡
3.添加付款方式（信用卡）
4. 添加初始积分（建议至少 5 美元）
5.可选择设置自动充值

**定价说明：**
- 模型具有不同的每个代币成本
- 请参阅 https://openrouter.ai/perplexity 了解困惑定价
- 监控 https://openrouter.ai/activity 的使用情况

### 第 3 步：生成 API 密钥

1. 前往https://openrouter.ai/keys
2. 点击“创建密钥”
3. 为您的密钥指定一个描述性名称（例如“perplexity-search-skill”）
4. 可选择设置安全使用限制
5.复制密钥（以`sk-or-v1-...`开头）
6. **重要**：安全保存此密钥 - 您无法再次查看它！

**安全提示：**
- 切勿公开分享您的 API 密钥
- 不要将密钥提交给版本控制
- 对不同的项目使用单独的密钥
- 设置使用限制以防止意外收费
- 定期轮换钥匙

### 步骤4：配置环境

您有两种设置 API 密钥的选项：

#### 选项 A：环境变量（推荐）

**Linux/macOS：**
```bash
export OPENROUTER_API_KEY='sk-or-v1-your-key-here'
```

要使其永久生效，请添加到您的 shell 配置文件中：
<<<代码块_1>>>

**Windows（PowerShell）：**
<<<代码块_2>>>

使其永久化：
<<<代码块_3>>>

#### 选项 B：.env 文件

在项目目录中创建一个 `.env` 文件：

<<<代码块_4>>>

或者使用安装脚本：
<<<代码块_5>>>

然后在运行脚本之前加载它：
<<<代码块_6>>>

**在脚本中使用 python-dotenv:**
```python
from dotenv import load_dotenv
load_dotenv()  # Loads .env file automatically

import os
api_key = os.environ.get("OPENROUTER_API_KEY")
```

### 步骤 5：安装依赖项

使用 uv 安装 LiteLLM：

```bash
uv pip install litellm
```

或者使用常规点：
```bash
pip install litellm
```

**可选依赖项：**
```bash
# For .env file support
uv pip install python-dotenv

# For additional features
uv pip install litellm[proxy]  # If using LiteLLM proxy server
```

### 第 6 步：验证设置

测试您的配置：

```bash
# Using the setup script
python scripts/setup_env.py --validate

# Or using the search script
python scripts/perplexity_search.py --check-setup
```

你应该看到：
```
✓ OPENROUTER_API_KEY is set (sk-or-v1-...xxxx)
✓ LiteLLM is installed (version X.X.X)
✓ Setup is complete! You're ready to use Perplexity Search.
```

### 第 7 步：测试您的第一次搜索

运行一个简单的测试查询：

```bash
python scripts/perplexity_search.py "What is CRISPR gene editing?"
```

预期输出：
```
================================================================================
ANSWER
================================================================================
CRISPR (Clustered Regularly Interspaced Short Palindromic Repeats) is a
revolutionary gene editing technology that allows precise modifications to DNA...
[detailed answer continues]
================================================================================
```

## 使用情况监控

### 检查您的使用情况

监控您的 OpenRouter 使用情况和成本：

1.访问https://openrouter.ai/activity
2.查看请求、令牌和成本
3. 按日期范围、型号或键进行过滤
4.导出使用数据进行分析

### 设置使用限制

防止意外收费：

1. 前往https://openrouter.ai/keys
2. 单击您的密钥
3. 设置“速率限制”（每分钟请求数）
4. 设置“支出限额”（最高总支出）
5. 如果需要，启用有限制的“自动充值”

**建议的开发限制：**
- 速率限制：每分钟 10-20 个请求
- 消费限额：10-50 美元，具体取决于使用情况

### 成本优化

降低成本的技巧：

1. **选择合适的型号**：使用Sonar进行简单查询，而不是Sonar Pro Search
2. **设置max_tokens**：使用`--max-tokens`参数限制响应长度
3. **批量查询**：尽可能组合多个简单问题
4. **监控使用情况**：在繁重的开发过程中每天检查成本
5. **使用缓存**：存储重复查询的结果

## 故障排除

### 错误：“未配置 OpenRouter API 密钥”

**原因**：环境变量未设置

**解决方案**：
```bash
# Check if variable is set
echo $OPENROUTER_API_KEY

# If empty, set it
export OPENROUTER_API_KEY='sk-or-v1-your-key-here'

# Or use setup script
python scripts/setup_env.py --api-key sk-or-v1-your-key-here
```

### 错误：“API 密钥无效”

**原因**：
- 密钥被删除或撤销
- 密钥已过期
- 键中的拼写错误
- 密钥格式错误

**解决方案**：
1. 验证 https://openrouter.ai/keys 处的密钥
2. 检查是否有多余的空格或引号
3. 如果需要，生成新密钥
4. 确保密钥以 `sk-or-v1-` 开头

### 错误：“积分不足”

**原因**：OpenRouter 帐户已用完积分

**解决方案**：
1. 前往https://openrouter.ai/account
2. 单击“制作人员”选项卡
3.添加更多学分
4.考虑启用自动充值

### 错误：“超出速率限制”

**原因**：短时间内请求过多

**解决方案**：
1. 等待几秒钟再重试
2. 提高 https://openrouter.ai/keys 的速率限制
3. 在代码中实现指数退避
4.批量请求或降低频率

### 错误：“找不到模型”

**原因**：型号名称不正确或型号不再可用

**解决方案**：
1. 在https://openrouter.ai/models查看可用型号
2. 使用正确的格式：`openrouter/perplexity/sonar-pro`
3. 验证模型是否仍受支持

### 错误：“未安装 LiteLLM”

**原因**：未安装 LiteLLM 软件包

**解决方案**：
```bash
uv pip install litellm
```

### LiteLLM 导入错误

**原因**：Python路径问题或版本冲突

**解决方案**：
1. 验证安装：`pip list | grep litellm`
2. 重新安装：`uv pip install --force-reinstall litellm`
3. 检查Python版本：`python --version`（需要3.8+）
4.使用虚拟环境避免冲突

## 高级配置

### 使用多个键

对于不同的项目或团队成员：

```bash
# Project 1
export OPENROUTER_API_KEY='sk-or-v1-project1-key'

# Project 2
export OPENROUTER_API_KEY='sk-or-v1-project2-key'
```

或者使用不同目录中的 .env 文件。

### 自定义基本 URL

如果使用 OpenRouter 代理或自定义端点：

```python
from litellm import completion

response = completion(
    model="openrouter/perplexity/sonar-pro",
    messages=[{"role": "user", "content": "query"}],
    api_base="https://custom-endpoint.com/v1"  # Custom URL
)
```

### 请求标头

添加自定义标头以进行跟踪：

```python
from litellm import completion

response = completion(
    model="openrouter/perplexity/sonar-pro",
    messages=[{"role": "user", "content": "query"}],
    extra_headers={
        "HTTP-Referer": "https://your-app.com",
        "X-Title": "Your App Name"
    }
)
```

### 超时配置

为长时间运行的查询设置自定义超时：

```python
from litellm import completion

response = completion(
    model="openrouter/perplexity/sonar-pro-search",
    messages=[{"role": "user", "content": "complex query"}],
    timeout=120  # 120 seconds timeout
)
```

## 安全最佳实践

### API密钥管理

1. **从不提交密钥**：将 `.env` 添加到 `.gitignore`
2. **使用密钥轮换**：每 3-6 个月轮换一次密钥
3. **单独的密钥**：用于开发/登台/生产的不同密钥
4. **监控使用情况**：检查是否有未经授权的访问
5. **设置限制**：配置支出和费率限制

### .gitignore 模板

添加到您的`.gitignore`：
```
# Environment variables
.env
.env.local
.env.*.local

# API keys
*api_key*
*apikey*
*.key

# Sensitive configs
config/secrets.yaml
```

### 密钥撤销

如果密钥被泄露：

1. 立即前往https://openrouter.ai/keys
2. 单击受损密钥上的“删除”
3. 生成新密钥
4. 使用旧密钥更新所有应用程序
5. 检查使用日志是否存在未经授权的访问
6. 如果需要，联系 OpenRouter 支持

## 常见问题解答

**问：通过 OpenRouter 使用 Perplexity 需要多少钱？**

答：价格因型号而异。 Sonar 最便宜（每次查询约 0.001-0.002 美元），Sonar Pro 中等（约 0.002-0.005 美元），Sonar Pro Search 最贵（每次查询约 0.02-0.05 美元+）。请参阅 https://openrouter.ai/perplexity 了解确切的定价。

**问：我需要单独的 Perplexity API 密钥吗？**

答：不！ OpenRouter 仅使用您的 OpenRouter 密钥即可访问 Perplexity 模型。

**问：除了 Perplexity 之外，我可以将 OpenRouter 用于其他模型吗？**

答：是的！ OpenRouter 可以通过相同的 API 密钥访问来自 OpenAI、Anthropic、Google、Meta 等的 100 多个模型。

**问：有免费套餐吗？**

答：OpenRouter 需要付费，但提供非常有竞争力的价格。最初的 5 美元信用额应该可以用于广泛的测试。

**问：如何取消我的 OpenRouter 帐户？**

答：联系 OpenRouter 支持。请注意，未使用的积分可能无法退还。

**问：我可以在生产应用程序中使用 OpenRouter 吗？**

答：是的，OpenRouter 专为生产用途而设计，具有强大的基础设施、SLA 和企业支持。

## 资源

**官方文档：**
- OpenRouter：https://openrouter.ai/docs
- 困惑模型：https://openrouter.ai/perplexity
- LiteLLM：https://docs.litellm.ai/

**账户管理：**
- 仪表板：https://openrouter.ai/account
- API 密钥：https://openrouter.ai/keys
- 用法：https://openrouter.ai/activity
- 计费：https://openrouter.ai/credits

**社区：**
- OpenRouter Discord：https://discord.gg/openrouter
- GitHub 问题：https://github.com/OpenRouter
- LiteLLM GitHub：https://github.com/BerriAI/litellm

## 总结

设置 OpenRouter 以进行 Perplexity 访问涉及：

1. 在 https://openrouter.ai 创建帐户
2. 添加付款方式和积分
3. 在 https://openrouter.ai/keys 生成 API 密钥
4. 设置`OPENROUTER_API_KEY`环境变量
5. 安装 LiteLLM：`uv pip install litellm`
6. 验证设置：`python scripts/setup_env.py --validate`
7. 开始搜索：`python scripts/perplexity_search.py "your query"`

定期监控使用情况和成本，以优化您的支出并确保安全。