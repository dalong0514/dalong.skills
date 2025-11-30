<!-- 此文件由机器翻译自 secrets.md -->

# 秘密和环境变量

## 创建秘密

### 通过仪表板

在 https://modal.com/secrets 创建机密

模板可用于：
- 数据库凭证（Postgres、MongoDB）
- 云提供商（AWS、GCP、Azure）
- ML 平台（权重和偏差、拥抱脸部）
- 还有更多

### 通过 CLI

```bash
# Create secret with key-value pairs
modal secret create my-secret KEY1=value1 KEY2=value2

# Use environment variables
modal secret create db-secret PGHOST=uri PGPASSWORD="$PGPASSWORD"

# List secrets
modal secret list

# Delete secret
modal secret delete my-secret
```

### 以编程方式

来自字典：

<<<代码块_1>>>

从 .env 文件：

<<<代码块_2>>>

## 使用秘密

将秘密注入到函数中：

<<<代码块_3>>>

### 多个秘密

<<<代码块_4>>>

如果密钥冲突，后面的秘密会覆盖前面的秘密。

## 环境变量

### 保留的运行时变量

**所有容器**：
- `MODAL_CLOUD_PROVIDER` - 云提供商 (AWS/GCP/OCI)
- `MODAL_IMAGE_ID` - 图像 ID
- `MODAL_REGION` - 区域标识符（例如 us-east-1）
- `MODAL_TASK_ID` - 容器任务 ID

**功能容器**：
- `MODAL_ENVIRONMENT` - 模态环境名称
- `MODAL_IS_REMOTE` - 在远程容器中设置为“1”
- `MODAL_IDENTITY_TOKEN` - 用于函数身份的 OIDC 令牌

**沙箱容器**：
- `MODAL_SANDBOX_ID` - 沙箱 ID

### 设置环境变量

通过图片：

<<<代码块_5>>>

通过秘密：

<<<代码块_6>>>

## 常见的秘密模式

### AWS 凭证

```python
aws_secret = modal.Secret.from_name("my-aws-secret")

@app.function(secrets=[aws_secret])
def use_aws():
    import boto3
    s3 = boto3.client('s3')
    # AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY automatically used
```

### 拥抱脸令牌

```python
hf_secret = modal.Secret.from_name("huggingface")

@app.function(secrets=[hf_secret])
def download_model():
    from transformers import AutoModel
    # HF_TOKEN automatically used for authentication
    model = AutoModel.from_pretrained("private-model")
```

### 数据库凭证

```python
db_secret = modal.Secret.from_name("postgres-creds")

@app.function(secrets=[db_secret])
def query_db():
    import psycopg2
    conn = psycopg2.connect(
        host=os.environ["PGHOST"],
        port=os.environ["PGPORT"],
        user=os.environ["PGUSER"],
        password=os.environ["PGPASSWORD"],
    )
```

## 最佳实践

1. **永远不要对秘密进行硬编码** - 始终使用模态秘密
2. **使用特定的秘密** - 为不同的目的创建单独的秘密
3. **定期轮换秘密** - 定期更新秘密
4. **最小范围** - 仅将秘密附加到需要它们的函数
5. **环境特定** - 对开发/登台/生产使用不同的秘密

## 安全说明

- 秘密在静态时被加密
- 仅适用于明确请求它们的函数
- 未记录或暴露在仪表板中
- 可以适用于特定环境