<!-- 此文件由机器翻译自 distributed_training.md -->

# 分布式训练 - 综合指南

## 概述

PyTorch Lightning 提供了多种跨多个 GPU、节点和机器有效训练大型模型的策略。根据模型大小和硬件配置选择正确的策略。

## 策略选择指南

### 何时使用每种策略

**定期培训（单台设备）**
- 模型大小：适合单个 GPU 内存的任何大小
- 使用案例：原型设计、小型模型、调试

**DDP（分布式数据并行）**
- 模型大小：<500M参数（例如ResNet50~80M参数）
- 何时：权重、激活、优化器状态和梯度都适合 GPU 内存
- 目标：跨多个 GPU 扩展批量大小和速度
- 最适合：最标准的深度学习模型

**FSDP（完全分片数据并行）**
- 模型大小：500M+参数（例如，BERT-Large、GPT 等大型变压器）
- 何时：模型不适合单个 GPU 内存
- 推荐对象：刚刚接触并行模型或从 DDP 迁移的用户
- 功能：激活检查点、CPU 参数卸载

**深速**
- 模型大小：500M+参数
- 何时：需要尖端功能或已经熟悉 DeepSpeed
- 特点：CPU/磁盘参数卸载、分布式检查点、细粒度控制
- 权衡：更复杂的配置

## DDP（分布式数据并行）

### 基本用法

```python
# Single GPU
trainer = L.Trainer(accelerator="gpu", devices=1)

# Multi-GPU on single node (automatic DDP)
trainer = L.Trainer(accelerator="gpu", devices=4)

# Explicit DDP strategy
trainer = L.Trainer(strategy="ddp", accelerator="gpu", devices=4)
```

### 多节点 DDP

<<<代码块_1>>>

### DDP 配置

<<<代码块_2>>>

### DDP 生成

当 `ddp` 导致问题时使用（速度较慢但更兼容）：

<<<代码块_3>>>

### DDP 最佳实践

1. **批量大小：** 乘以 GPU 数量
   <<<代码块_4>>>

2. **学习率：** 通常与批量大小成比例
   <<<代码块_5>>>

3. **同步：** 使用 `sync_dist=True` 作为指标
   <<<代码块_6>>>

4. **特定于等级的操作：** 仅在主流程中使用装饰器
   ```python
   from lightning.pytorch.utilities import rank_zero_only

   @rank_zero_only
   def save_results(self):
       # Only runs on main process (rank 0)
       torch.save(self.results, "results.pt")
   ```

## FSDP（完全分片数据并行）

### 基本用法

```python
trainer = L.Trainer(
    strategy="fsdp",
    accelerator="gpu",
    devices=4
)
```

### FSDP 配置

```python
from lightning.pytorch.strategies import FSDPStrategy
import torch.nn as nn

trainer = L.Trainer(
    strategy=FSDPStrategy(
        # Sharding strategy
        sharding_strategy="FULL_SHARD",  # or "SHARD_GRAD_OP", "NO_SHARD", "HYBRID_SHARD"

        # Activation checkpointing (save memory)
        activation_checkpointing_policy={nn.TransformerEncoderLayer},

        # CPU offloading (save GPU memory, slower)
        cpu_offload=False,

        # Mixed precision
        mixed_precision=True,

        # Wrap policy (auto-wrap layers)
        auto_wrap_policy=None
    ),
    accelerator="gpu",
    devices=8,
    precision="bf16-mixed"
)
```

### 分片策略

**FULL_SHARD（默认）**
- 分片优化器状态、梯度和参数
- 最大程度节省内存
- 更多的通信开销

**SHARD_GRAD_OP**
- 仅分片优化器状态和梯度
- 参数保存在所有设备上
- 节省的内存更少，但速度更快

**NO_SHARD**
- 无分片（相当于DDP）
- 用于比较或不需要分片时

**混合_碎片**
- 结合节点内的 FULL_SHARD 和跨节点的 NO_SHARD
- 适合多节点设置

### 激活检查点

内存的交易计算：

```python
from lightning.pytorch.strategies import FSDPStrategy
import torch.nn as nn

# Checkpoint specific layer types
trainer = L.Trainer(
    strategy=FSDPStrategy(
        activation_checkpointing_policy={
            nn.TransformerEncoderLayer,
            nn.TransformerDecoderLayer
        }
    )
)
```

### CPU 卸载

不使用时将参数卸载到 CPU：

```python
trainer = L.Trainer(
    strategy=FSDPStrategy(
        cpu_offload=True  # Slower but saves GPU memory
    ),
    accelerator="gpu",
    devices=4
)
```

### 带有大型模型的 FSDP

```python
from lightning.pytorch.strategies import FSDPStrategy
import torch.nn as nn

class LargeTransformer(L.LightningModule):
    def __init__(self):
        super().__init__()
        self.transformer = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d_model=4096, nhead=32),
            num_layers=48
        )

    def configure_sharded_model(self):
        # Called by FSDP to wrap model
        pass

# Train
trainer = L.Trainer(
    strategy=FSDPStrategy(
        activation_checkpointing_policy={nn.TransformerEncoderLayer},
        cpu_offload=False,
        sharding_strategy="FULL_SHARD"
    ),
    accelerator="gpu",
    devices=8,
    precision="bf16-mixed",
    max_epochs=10
)

model = LargeTransformer()
trainer.fit(model, datamodule=dm)
```

## 深速

### 安装

```bash
pip install deepspeed
```

### 基本用法

```python
trainer = L.Trainer(
    strategy="deepspeed_stage_2",  # or "deepspeed_stage_3"
    accelerator="gpu",
    devices=4,
    precision="16-mixed"
)
```

### DeepSpeed 阶段

**阶段 1：优化器状态分片**
- 分片优化器状态
- 适度节省内存

```python
trainer = L.Trainer(strategy="deepspeed_stage_1")
```

**第二阶段：优化器+梯度分片**
- 分片优化器状态和梯度
- 良好的内存节省

```python
trainer = L.Trainer(strategy="deepspeed_stage_2")
```

**第 3 阶段：完整模型分片 (ZeRO-3)**
- 分片优化器状态、梯度和模型参数
- 最大程度节省内存
- 可以训练非常大的模型

```python
trainer = L.Trainer(strategy="deepspeed_stage_3")
```

**第二阶段卸载**
- 卸载到 CPU 或 NVMe

```python
trainer = L.Trainer(strategy="deepspeed_stage_2_offload")
trainer = L.Trainer(strategy="deepspeed_stage_3_offload")
```

### DeepSpeed 配置文件

对于细粒度控制：

```python
from lightning.pytorch.strategies import DeepSpeedStrategy

# Create config file: ds_config.json
config = {
    "zero_optimization": {
        "stage": 3,
        "offload_optimizer": {
            "device": "cpu",
            "pin_memory": True
        },
        "offload_param": {
            "device": "cpu",
            "pin_memory": True
        },
        "overlap_comm": True,
        "contiguous_gradients": True,
        "sub_group_size": 1e9,
        "reduce_bucket_size": "auto",
        "stage3_prefetch_bucket_size": "auto",
        "stage3_param_persistence_threshold": "auto",
        "stage3_max_live_parameters": 1e9,
        "stage3_max_reuse_distance": 1e9
    },
    "fp16": {
        "enabled": True,
        "loss_scale": 0,
        "initial_scale_power": 16,
        "loss_scale_window": 1000,
        "hysteresis": 2,
        "min_loss_scale": 1
    },
    "gradient_clipping": 1.0,
    "train_batch_size": "auto",
    "train_micro_batch_size_per_gpu": "auto"
}

trainer = L.Trainer(
    strategy=DeepSpeedStrategy(config=config),
    accelerator="gpu",
    devices=8,
    precision="16-mixed"
)
```

### DeepSpeed 最佳实践

1. **对 <10B 参数的模型使用第 2 阶段**
2. **对于 >10B 参数的模型使用第 3 阶段**
3. **GPU内存不足时启用卸载**
4. **调整 `reduce_bucket_size` 以提高通信效率**

## 比较表

|特色 |完税后交货 | FSDP |深速|
|--------|-----|------|------------|
|型号尺寸| <500M 参数 | 500M+ 参数 | 500M+ 参数 |
|内存效率 |低|高|非常高 |
|速度|最快|快|快|
|设置复杂性 |简单|中等|复杂|
|卸载 |没有 |中央处理器| CPU + 磁盘 |
|最适合 |标准型号|大型号|非常大的模型|
|配置|最小 |中等|广泛 |

## 混合精度训练

使用混合精度来加速训练并节省内存：

```python
# FP16 mixed precision
trainer = L.Trainer(precision="16-mixed")

# BFloat16 mixed precision (A100, H100)
trainer = L.Trainer(precision="bf16-mixed")

# Full precision (default)
trainer = L.Trainer(precision="32-true")

# Double precision
trainer = L.Trainer(precision="64-true")
```

### 不同策略的混合精度

```python
# DDP + FP16
trainer = L.Trainer(
    strategy="ddp",
    accelerator="gpu",
    devices=4,
    precision="16-mixed"
)

# FSDP + BFloat16
trainer = L.Trainer(
    strategy="fsdp",
    accelerator="gpu",
    devices=8,
    precision="bf16-mixed"
)

# DeepSpeed + FP16
trainer = L.Trainer(
    strategy="deepspeed_stage_2",
    accelerator="gpu",
    devices=4,
    precision="16-mixed"
)
```

## 多节点训练

### 胡言乱语

```bash
#!/bin/bash
#SBATCH --nodes=4
#SBATCH --gpus-per-node=4
#SBATCH --time=24:00:00

srun python train.py
```

```python
# train.py
trainer = L.Trainer(
    strategy="ddp",
    accelerator="gpu",
    devices=4,
    num_nodes=4
)
```

### 手动多节点设置

节点 0（主节点）：
```bash
python train.py --num_nodes=2 --node_rank=0 --master_addr=192.168.1.1 --master_port=12345
```

节点1：
```bash
python train.py --num_nodes=2 --node_rank=1 --master_addr=192.168.1.1 --master_port=12345
```

```python
# train.py
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--num_nodes", type=int, default=1)
parser.add_argument("--node_rank", type=int, default=0)
parser.add_argument("--master_addr", type=str, default="localhost")
parser.add_argument("--master_port", type=int, default=12345)
args = parser.parse_args()

trainer = L.Trainer(
    strategy="ddp",
    accelerator="gpu",
    devices=4,
    num_nodes=args.num_nodes
)
```

## 常见模式
### 使用 DDP 进行梯度累积

```python
# Simulate larger batch size
trainer = L.Trainer(
    strategy="ddp",
    accelerator="gpu",
    devices=4,
    accumulate_grad_batches=4  # Effective batch size = batch_size * devices * 4
)
```

### 具有分布式训练的模型检查点

```python
from lightning.pytorch.callbacks import ModelCheckpoint

checkpoint_callback = ModelCheckpoint(
    monitor="val_loss",
    save_top_k=3,
    mode="min"
)

trainer = L.Trainer(
    strategy="ddp",
    accelerator="gpu",
    devices=4,
    callbacks=[checkpoint_callback]
)
```

### 分布式训练的可重复性

```python
import lightning as L

L.seed_everything(42, workers=True)

trainer = L.Trainer(
    strategy="ddp",
    accelerator="gpu",
    devices=4,
    deterministic=True
)
```

## 故障排除

### NCCL 超时

增加慢速网络的超时：

```python
import os
os.environ["NCCL_TIMEOUT"] = "3600"  # 1 hour

trainer = L.Trainer(strategy="ddp", accelerator="gpu", devices=4)
```

### CUDA 内存不足

解决方案：
1.启用梯度检查点
2. 减少批量
3.使用FSDP或DeepSpeed
4.启用CPU卸载
5.使用混合精度

```python
# Option 1: Gradient checkpointing
class MyModel(L.LightningModule):
    def __init__(self):
        super().__init__()
        self.model = MyTransformer()
        self.model.gradient_checkpointing_enable()

# Option 2: Smaller batch size
dm = MyDataModule(batch_size=16)  # Reduce from 32

# Option 3: FSDP with offloading
trainer = L.Trainer(
    strategy=FSDPStrategy(cpu_offload=True),
    precision="bf16-mixed"
)

# Option 4: Gradient accumulation
trainer = L.Trainer(accumulate_grad_batches=4)
```

### 分布式采样器问题

Lightning 自动处理 DistributedSampler：

```python
# Don't do this
from torch.utils.data import DistributedSampler
sampler = DistributedSampler(dataset)  # Lightning does this automatically

# Just use shuffle
train_loader = DataLoader(dataset, batch_size=32, shuffle=True)
```

### 通信开销

减少与较大的 `find_unused_parameters` 的通信：

```python
trainer = L.Trainer(
    strategy=DDPStrategy(find_unused_parameters=False),
    accelerator="gpu",
    devices=4
)
```

## 最佳实践

### 1. 从单 GPU 开始
在缩放之前在单个 GPU 上测试您的代码：

```python
# Debug on single GPU
trainer = L.Trainer(accelerator="gpu", devices=1, fast_dev_run=True)

# Then scale to multiple GPUs
trainer = L.Trainer(accelerator="gpu", devices=4, strategy="ddp")
```

### 2. 使用适当的策略
- <500M 参数：使用 DDP
- 500M-10B参数：使用FSDP
- >10B 参数：使用 DeepSpeed Stage 3

### 3.启用混合精度
始终对现代 GPU 使用混合精度：

```python
trainer = L.Trainer(precision="bf16-mixed")  # A100, H100
trainer = L.Trainer(precision="16-mixed")    # V100, T4
```

### 4. 缩放超参数
缩放时调整学习率和批量大小：

```python
# Linear scaling rule
lr = base_lr * num_gpus
```

### 5. 同步指标
始终同步分布式训练中的指标：

```python
self.log("val_loss", loss, sync_dist=True)
```

### 6. 使用零阶运算
仅在主进程上进行文件 I/O 和昂贵的操作：

```python
from lightning.pytorch.utilities import rank_zero_only

@rank_zero_only
def save_predictions(self):
    torch.save(self.predictions, "predictions.pt")
```

### 7. 定期检查
保存检查点以从故障中恢复：

```python
checkpoint_callback = ModelCheckpoint(
    save_top_k=3,
    save_last=True,  # Always save last for resuming
    every_n_epochs=5
)
```