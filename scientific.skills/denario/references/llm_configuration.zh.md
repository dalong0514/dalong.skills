<!-- 此文件由机器翻译自 llm_configuration.md -->

#LLM API配置

## 概述

Denario 需要来自受支持的 LLM 提供商的 API 凭证来为其多代理研究系统提供支持。该系统基于AG2和LangGraph构建，支持多个LLM后端。

## 支持的法学硕士提供商

### 谷歌顶点人工智能
- 与 Google Vertex AI 平台完全集成
- 支持Gemini和PaLM模型
- 需要 Google Cloud 项目设置

### 开放人工智能
- GPT-4、GPT-3.5 和其他 OpenAI 模型
- 直接API集成

### 其他提供商
- 任何与 AG2/LangGraph 框架兼容的法学硕士
- 人类克劳德（通过兼容接口）
- Azure OpenAI
- 自定义模型端点

## 获取API密钥

### 谷歌顶点人工智能

1. **创建Google云项目**
   - 导航到 [Google Cloud Console](https://console.cloud.google.com/)
   - 创建一个新项目或选择现有项目

2. **启用 Vertex AI API**
   - 转到“API 和服务”→“库”
   - 搜索“Vertex AI API”
   - 单击“启用”

3. **创建服务帐户**
   - 导航到“IAM 和管理”→“服务帐户”
   - 创建具有 Vertex AI 权限的服务帐户
   - 下载 JSON 密钥文件

4. **设置身份验证**
   ```bash
   export GOOGLE_APPLICATION_CREDENTIALS="/path/to/service-account-key.json"
   ```

### 开放人工智能

1. **创建OpenAI账户**
   - 访问 [platform.openai.com](https://platform.openai.com/)
   - 注册或登录

2. **生成API密钥**
   - 导航至 API 密钥部分
   - 单击“创建新密钥”
   - 安全地复制和存储

3. **设置环境变量**
   <<<代码块_1>>>

## 存储 API 密钥

### 方法一：环境变量（推荐）

**Linux/macOS：**
<<<代码块_2>>>

添加到 `~/.bashrc`、`~/.zshrc` 或 `~/.bash_profile` 以实现持久性。

**Windows：**
<<<代码块_3>>>

或者使用系统属性→环境变量进行持久化。

### 方法 2：.env 文件

在项目目录中创建一个 `.env` 文件：

<<<代码块_4>>>

在Python中加载环境文件：

<<<代码块_5>>>

### 方法3：Docker环境文件

对于 Docker 部署，传递环境变量：

<<<代码块_6>>>

## Vertex AI 详细设置

### 先决条件
- 已启用结算功能的 Google Cloud 帐户
- 安装 gcloud CLI（可选但推荐）

### 逐步配置

1. **安装 Google Cloud SDK（如果不使用 Docker）**
   ```bash
   # Linux/macOS
   curl https://sdk.cloud.google.com | bash
   exec -l $SHELL
   gcloud init
   ```

2. **验证 gcloud**
   ```bash
   gcloud auth application-default login
   ```

3. **设置项目**
   ```bash
   gcloud config set project YOUR_PROJECT_ID
   ```

4. **启用所需的API**
   ```bash
   gcloud services enable aiplatform.googleapis.com
   gcloud services enable compute.googleapis.com
   ```

5. **创建服务帐户（替代 gcloud auth）**
   ```bash
   gcloud iam service-accounts create denario-service-account \
     --display-name="Denario AI Service Account"

   gcloud projects add-iam-policy-binding YOUR_PROJECT_ID \
     --member="serviceAccount:denario-service-account@YOUR_PROJECT_ID.iam.gserviceaccount.com" \
     --role="roles/aiplatform.user"

   gcloud iam service-accounts keys create credentials.json \
     --iam-account=denario-service-account@YOUR_PROJECT_ID.iam.gserviceaccount.com
   ```

6. **配置 denario 以使用 Vertex AI**
   ```python
   import os
   os.environ['GOOGLE_CLOUD_PROJECT'] = 'YOUR_PROJECT_ID'
   os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = '/path/to/credentials.json'

   from denario import Denario
   den = Denario(project_dir="./research")
   ```

## 型号选择

配置 denario 用于不同任务的模型：

```python
# In your code
from denario import Denario

# Example configuration (if supported by denario API)
den = Denario(
    project_dir="./project",
    # Model configuration may vary based on denario version
)
```

检查 denario 的文档以了解特定模型选择 API。

## 成本管理

### 监控成本

- **OpenAI**：在 [platform.openai.com/usage](https://platform.openai.com/usage) 跟踪使用情况
- **Google Cloud**：在 Cloud Console 中监控 → 结算
- 设置账单提醒以避免意外收费

### 成本优化技巧

1. **使用适当的模型层**
   - GPT-3.5 用于更简单的任务
   - 用于复杂推理的 GPT-4

2. **批量操作**
   - 在单个会话中处理多项研究任务

3. **缓存结果**
   - 尽可能重用生成的想法、方法和结果

4. **设置代币限制**
   - 配置最大令牌使用量以控制成本

## 安全最佳实践

### 不要将 API 密钥提交给版本控制

添加到 `.gitignore`：
```gitignore
.env
*.json  # If storing credentials
credentials.json
service-account-key.json
```

### 定期轮换钥匙
- 定期生成新的 API 密钥
- 轮换后撤销旧密钥

### 使用最小权限访问
- 仅向服务帐户授予必要的权限
- 使用单独的密钥进行开发和生产

### 加密敏感文件
- 将凭证文件存储在加密卷中
- 使用云秘密管理服务进行生产

## 故障排除

### “未找到 API 密钥”错误
- 验证环境变量是否已设置：`echo $OPENAI_API_KEY`
- 检查`.env`文件是否位于正确的目录中
- 确保在导入 denario 之前调用 `load_dotenv()`

### Vertex AI 身份验证失败
- 验证 `GOOGLE_APPLICATION_CREDENTIALS` 指向有效的 JSON 文件
- 检查服务帐户是否具有所需的权限
- 确保在 Google Cloud 项目中启用 API

### 速率限制问题
- 实施指数退避
- 减少并发请求
- 如果需要升级API计划

### Docker环境变量问题
- 使用`docker run --env-file .env`传递环境
- 使用 `-v` 标志挂载凭证文件
- 检查容器内的环境：`docker exec <container> env`