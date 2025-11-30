<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：莫代尔
描述：使用无服务器容器、GPU 和自动缩放在云中运行 Python 代码。在部署 ML 模型、运行批处理作业、调度计算密集型任务或提供需要 GPU 加速或动态扩展的 API 时使用。
---

# 模态

## 概述

Modal 是一个无服务器平台，用于以最少的配置在云中运行 Python 代码。在强大的 GPU 上执行功能，自动扩展到数千个容器，并且只需为使用的计算付费。

Modal 特别适合 AI/ML 工作负载、高性能批处理、计划作业、GPU 推理和无服务器 API。通过 https://modal.com 免费注册并获得每月 30 美元的积分。

## 何时使用此技能

使用模态：
- 部署和服务 ML 模型（LLM、图像生成、嵌入模型）
- 运行 GPU 加速计算（训练、推理、渲染）
- 并行批处理大型数据集
- 调度计算密集型作业（日常数据处理、模型训练）
- 构建需要自动扩展的无服务器API
- 需要分布式计算或专用硬件的科学计算

## 身份验证和设置

Modal 需要通过 API 令牌进行身份验证。

### 初始设置

```bash
# Install Modal
uv uv pip install modal

# Authenticate (opens browser for login)
modal token new
```

这将创建一个存储在 `~/.modal.toml` 中的令牌。令牌对所有 Modal 操作进行身份验证。

### 验证设置

<<<代码块_1>>>

运行方式：`modal run script.py`

## 核心能力

Modal 通过在容器中运行的函数提供无服务器 Python 执行。以声明方式定义计算要求、依赖关系和扩展行为。

### 1. 定义容器镜像

使用模态图像指定函数的依赖关系和环境。

<<<代码块_2>>>

**常见模式：**
- 安装Python包：`.uv_pip_install("pandas", "scikit-learn")`
- 安装系统包：`.apt_install("ffmpeg", "git")`
- 使用现有的 Docker 镜像：`modal.Image.from_registry("nvidia/cuda:12.1.0-base")`
- 添加本地代码：`.add_local_python_source("my_module")`

有关全面的映像构建文档，请参阅 `references/images.md`。

### 2. 创建函数

使用 `@app.function()` 装饰器定义在云中运行的函数。

<<<代码块_3>>>

**调用功能：**
<<<代码块_4>>>

运行：`modal run script.py`

请参阅 `references/functions.md` 了解函数模式、部署和参数处理。

### 3. 请求 GPU

将 GPU 连接到函数以加速计算。

<<<代码块_5>>>

**可用的 GPU 类型：**
- `T4`、`L4` - 经济高效的推理
- `A10`、`A100`、`A100-80GB` - 标准训练/推理
- `L40S` - 卓越的成本/性能平衡 (48GB)
- `H100`、`H200` - 高性能训练
- `B200` - 旗舰性能（最强大）

**请求多个 GPU：**
<<<代码块_6>>>

请参阅 `references/gpu.md` 了解 GPU 选择指南、CUDA 设置和多 GPU 配置。

### 4.配置资源

请求 CPU 核心、内存和磁盘的功能。

```python
@app.function(
    cpu=8.0,           # 8 physical cores
    memory=32768,      # 32 GiB RAM
    ephemeral_disk=10240  # 10 GiB disk
)
def memory_intensive_task():
    pass
```

默认分配：0.125 个 CPU 核心，128 MiB 内存。根据预订或实际使用情况（以较高者为准）计费。

请参阅 `references/resources.md` 了解资源限制和计费详细信息。

### 5. 自动扩展

Modal 可根据需求自动扩展功能，从零到数千个容器。

**并行处理输入：**
```python
@app.function()
def analyze_sample(sample_id: int):
    # Process single sample
    return result

@app.local_entrypoint()
def main():
    sample_ids = range(1000)
    # Automatically parallelized across containers
    results = list(analyze_sample.map(sample_ids))
```

**配置自动缩放：**
```python
@app.function(
    max_containers=100,      # Upper limit
    min_containers=2,        # Keep warm
    buffer_containers=5      # Idle buffer for bursts
)
def inference():
    pass
```

有关自动缩放配置、并发和缩放限制，请参阅 `references/scaling.md`。

### 6. 持久存储数据

使用卷在函数调用之间进行持久存储。

```python
volume = modal.Volume.from_name("my-data", create_if_missing=True)

@app.function(volumes={"/data": volume})
def save_results(data):
    with open("/data/results.txt", "w") as f:
        f.write(data)
    volume.commit()  # Persist changes
```

卷在运行之间保留数据、存储模型权重、缓存数据集以及在函数之间共享数据。

有关卷管理、提交和缓存模式，请参阅 `references/volumes.md`。

### 7. 管理秘密

使用模态机密安全地存储 API 密钥和凭据。

```python
@app.function(secrets=[modal.Secret.from_name("huggingface")])
def download_model():
    import os
    token = os.environ["HF_TOKEN"]
    # Use token for authentication
```

**在 Modal 仪表板或通过 CLI 创建机密：**
```bash
modal secret create my-secret KEY=value API_TOKEN=xyz
```

有关秘密管理和身份验证模式，请参阅 `references/secrets.md`。

### 8. 部署 Web 端点

使用 `@modal.web_endpoint()` 提供 HTTP 端点、API 和 Webhook。

```python
@app.function()
@modal.web_endpoint(method="POST")
def predict(data: dict):
    # Process request
    result = model.predict(data["input"])
    return {"prediction": result}
```

**部署：**
```bash
modal deploy script.py
```

Modal 为端点提供 HTTPS URL。

请参阅 `references/web-endpoints.md` 了解 FastAPI 集成、流式传输、身份验证和 WebSocket 支持。

### 9. 安排工作

使用 cron 表达式按计划运行函数。

```python
@app.function(schedule=modal.Cron("0 2 * * *"))  # Daily at 2 AM
def daily_backup():
    # Backup data
    pass

@app.function(schedule=modal.Period(hours=4))  # Every 4 hours
def refresh_cache():
    # Update cache
    pass
```
计划的功能自动运行，无需手动调用。

请参阅 `references/scheduled-jobs.md` 了解 cron 语法、时区配置和监控。

## 常见工作流程

### 部署 ML 模型进行推理

```python
import modal

# Define dependencies
image = modal.Image.debian_slim().uv_pip_install("torch", "transformers")
app = modal.App("llm-inference", image=image)

# Download model at build time
@app.function()
def download_model():
    from transformers import AutoModel
    AutoModel.from_pretrained("bert-base-uncased")

# Serve model
@app.cls(gpu="L40S")
class Model:
    @modal.enter()
    def load_model(self):
        from transformers import pipeline
        self.pipe = pipeline("text-classification", device="cuda")

    @modal.method()
    def predict(self, text: str):
        return self.pipe(text)

@app.local_entrypoint()
def main():
    model = Model()
    result = model.predict.remote("Modal is great!")
    print(result)
```

### 批量处理大型数据集

```python
@app.function(cpu=2.0, memory=4096)
def process_file(file_path: str):
    import pandas as pd
    df = pd.read_csv(file_path)
    # Process data
    return df.shape[0]

@app.local_entrypoint()
def main():
    files = ["file1.csv", "file2.csv", ...]  # 1000s of files
    # Automatically parallelized across containers
    for count in process_file.map(files):
        print(f"Processed {count} rows")
```

### 在 GPU 上训练模型

```python
@app.function(
    gpu="A100:2",      # 2x A100 GPUs
    timeout=3600       # 1 hour timeout
)
def train_model(config: dict):
    import torch
    # Multi-GPU training code
    model = create_model(config)
    train(model)
    return metrics
```

## 参考文档

特定功能的详细文档：

- **`references/getting-started.md`** - 身份验证、设置、基本概念
- **`references/images.md`** - 镜像构建、依赖项、Dockerfile
- **`references/functions.md`** - 函数模式、部署、参数
- **`references/gpu.md`** - GPU 类型、CUDA、多 GPU 配置
- **`references/resources.md`** - CPU、内存、磁盘管理
- **`references/scaling.md`** - 自动缩放、并行执行、并发
- **`references/volumes.md`** - 持久存储、数据管理
- **`references/secrets.md`** - 环境变量、身份验证
- **`references/web-endpoints.md`** - API、webhooks、端点
- **`references/scheduled-jobs.md`** - Cron 作业，定期任务
- **`references/examples.md`** - 科学计算的常见模式

## 最佳实践

1. **在 `.uv_pip_install()` 中固定依赖项**，以实现可重现的构建
2. **使用适当的 GPU 类型** - L40S 用于推理，H100/A100 用于训练
3. **利用缓存** - 使用卷作为模型权重和数据集
4. **配置自动缩放** - 根据工作负载设置 `max_containers` 和 `min_containers`
5. **如果本地没有，则在函数体中导入包**
6. **使用 `.map()` 进行并行处理**而不是顺序循环
7. **安全地存储机密** - 切勿对 API 密钥进行硬编码
8. **监控成本** - 检查模态仪表板的使用情况和计费

## 故障排除

**“找不到模块”错误：**
- 使用 `.uv_pip_install("package-name")` 将包添加到映像
- 如果本地不可用，则在函数体内导入包

**未检测到 GPU：**
- 验证 GPU 规格：`@app.function(gpu="A100")`
- 检查 CUDA 可用性：`torch.cuda.is_available()`

**函数超时：**
- 增加超时：`@app.function(timeout=3600)`
- 默认超时为 5 分钟

**成交量变化不持久：**
- 写入文件后调用`volume.commit()`
- 验证函数装饰器中正确安装的卷

如需其他帮助，请参阅 https://modal.com/docs 处的 Modal 文档或加入 Modal Slack 社区。