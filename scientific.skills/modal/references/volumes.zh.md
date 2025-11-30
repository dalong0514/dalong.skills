<!-- 此文件由机器翻译自 volumes.md -->

# 模态体积

## 概述

Modal Volumes 为 Modal 应用程序提供高性能分布式文件系统。专为一次写入、多次读取的工作负载而设计，例如 ML 模型权重和分布式数据处理。

## 创建卷

### 通过 CLI

```bash
modal volume create my-volume
```

对于卷 v2（测试版）：
<<<代码块_1>>>

### 来自代码

<<<代码块_2>>>

## 使用体积

通过挂载点附加到函数：

<<<代码块_3>>>

## 提交和重新加载

### 提交

保留对音量的更改：

<<<代码块_4>>>

**后台提交**：Modal 每隔几秒以及在容器关闭时自动提交卷更改。

### 重新加载

从其他容器获取最新更改：

<<<代码块_5>>>

创建容器时，会安装最新的卷状态。需要重新加载才能查看其他容器的后续提交。

## 上传文件

### 批量上传（高效）

<<<代码块_6>>>

### 通过图像

```python
image = modal.Image.debian_slim().add_local_dir(
    local_path="/home/user/my_dir",
    remote_path="/app"
)

@app.function(image=image)
def process():
    # Files available at /app
    ...
```

## 下载文件

### 通过 CLI

```bash
modal volume get my-volume remote.txt local.txt
```

通过 CLI 的最大文件大小：无限制
通过仪表板的最大文件大小：16 MB

### 通过Python SDK

```python
vol = modal.Volume.from_name("my-volume")

for data in vol.read_file("path.txt"):
    print(data)
```

## 体积性能

### 卷 v1

最适合：
- <50,000 个文件（推荐）
- <500,000 个文件（硬限制）
- 顺序访问模式
- <5 个并发作者

### 卷 v2（测试版）

改进用于：
- 无限文件
- 数百名并发作家
- 随机访问模式
- 大文件（最大 1 TiB）

当前 v2 限制：
- 最大文件大小：1 TiB
- 每个目录的最大文件数：32,768
- 无限的目录深度

## 模型存储

### 保存模型权重

```python
volume = modal.Volume.from_name("model-weights", create_if_missing=True)
MODEL_DIR = "/models"

@app.function(volumes={MODEL_DIR: volume})
def train():
    model = train_model()
    save_model(f"{MODEL_DIR}/my_model.pt", model)
    volume.commit()
```

### 加载模型权重

```python
@app.function(volumes={MODEL_DIR: volume})
def inference(model_id: str):
    try:
        model = load_model(f"{MODEL_DIR}/{model_id}")
    except NotFound:
        volume.reload()  # Fetch latest models
        model = load_model(f"{MODEL_DIR}/{model_id}")
    return model.run(request)
```

## 模型检查点

在长时间训练作业期间保存检查点：

```python
volume = modal.Volume.from_name("checkpoints")
VOL_PATH = "/vol"

@app.function(
    gpu="A10G",
    timeout=2*60*60,  # 2 hours
    volumes={VOL_PATH: volume}
)
def finetune():
    from transformers import Seq2SeqTrainer, Seq2SeqTrainingArguments

    training_args = Seq2SeqTrainingArguments(
        output_dir=str(VOL_PATH / "model"),  # Checkpoints saved to Volume
        save_steps=100,
        # ... more args
    )

    trainer = Seq2SeqTrainer(model=model, args=training_args, ...)
    trainer.train()
```

即使训练中断，后台提交也可确保检查点持续存在。

## CLI 命令

```bash
# List files
modal volume ls my-volume

# Upload
modal volume put my-volume local.txt remote.txt

# Download
modal volume get my-volume remote.txt local.txt

# Copy within Volume
modal volume cp my-volume src.txt dst.txt

# Delete
modal volume rm my-volume file.txt

# List all volumes
modal volume list

# Delete volume
modal volume delete my-volume
```

## 临时卷

创建被垃圾收集的临时卷：

```python
with modal.Volume.ephemeral() as vol:
    sb = modal.Sandbox.create(
        volumes={"/cache": vol},
        app=my_app,
    )
    # Use volume
    # Automatically cleaned up when context exits
```

## 并发访问

### 并发读取

多个容器可以同时读取而不会出现问题。

### 并发写入

支持但是：
- 避免同时修改相同的文件
- 最后一次写入获胜（可能会丢失数据）
- v1：限制约 5 个并发写入者
- v2：支持数百个并发写入者

## 音量错误

###“音量忙”

文件打开时无法重新加载：

```python
# WRONG
f = open("/vol/data.txt", "r")
volume.reload()  # ERROR: volume busy
```

```python
# CORRECT
with open("/vol/data.txt", "r") as f:
    data = f.read()
# File closed before reload
volume.reload()
```

###“找不到文件”

请记住使用挂载点：

```python
# WRONG - file saved to local disk
with open("/xyz.txt", "w") as f:
    f.write("data")

# CORRECT - file saved to Volume
with open("/data/xyz.txt", "w") as f:
    f.write("data")
```

## 从 v1 升级到 v2

目前没有自动迁移。手动步骤：

1.创建新的v2卷
2. 使用 `cp` 或 `rsync` 复制数据
3.更新应用程序以使用新卷

```bash
modal volume create --version=2 my-volume-v2
modal shell --volume my-volume --volume my-volume-v2

# In shell:
cp -rp /mnt/my-volume/. /mnt/my-volume-v2/.
sync /mnt/my-volume-v2
```

警告：已部署的应用程序按 ID 引用卷。创建新Volume后重新部署。