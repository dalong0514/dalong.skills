<!-- 此文件由机器翻译自 gpu.md -->

# 模态上的 GPU 加速

## 快速入门

使用 `gpu` 参数在 GPU 上运行函数：

```python
import modal

image = modal.Image.debian_slim().pip_install("torch")
app = modal.App(image=image)

@app.function(gpu="A100")
def run():
    import torch
    assert torch.cuda.is_available()
```

## 可用的 GPU 类型

Modal 支持以下 GPU：

- `T4` - 入门级 GPU
- `L4` - 平衡性能和成本
- `A10` - 最多 4 个 GPU，总共 96 GB
- `A100` - 40GB 或 80GB 变体
- `A100-40GB` - 特定 40GB 变体
- `A100-80GB` - 特定 80GB 变体
- `L40S` - 48 GB，非常适合推理
- `H100` / `H100!` - 顶级 Hopper 架构
- `H200` - 改进了漏斗，具有更多内存
- `B200` - 最新的 Blackwell 架构

请参阅 https://modal.com/pricing 了解定价。

## GPU 数量

使用 `:n` 语法为每个容器请求多个 GPU：

<<<代码块_1>>>

支持的计数：
- B200、H200、H100、A100、L4、T4、L40S：最多 8 个 GPU（最多 1,536 GB）
- A10：最多 4 个 GPU（最高 96 GB）

注意：请求 >2 个 GPU 可能会导致更长的等待时间。

## GPU 选型指南

**用于推理（推荐）**：从 L40S 开始
- 卓越的性价比
- 48 GB 内存
- 适合 LLaMA、稳定扩散等。

**对于培训**：考虑 H100 或 A100
- 高计算吞吐量
- 大内存用于批处理

**对于内存密集型任务**：H200 或 A100-80GB
- 更多内存容量
- 更适合大型模型

## B200 GPU

NVIDIA 旗舰 Blackwell 芯片：

<<<代码块_2>>>

## H200 和 H100 GPU

Hopper 架构 GPU 具有出色的软件支持：

<<<代码块_3>>>

### H200 自动升级

Modal 可以将 `gpu="H100"` 升级到 H200，无需额外费用。 H200提供：
- 141 GB 内存（H100 为 80 GB）
- 4.8 TB/秒带宽（对比 3.35 TB/秒）

为了避免自动升级（例如，用于基准测试）：
<<<代码块_4>>>

## A100 GPU

Ampere 架构具有 40GB 或 80GB 型号：

<<<代码块_5>>>

## GPU 回退

指定具有后备功能的多种 GPU 类型：

<<<代码块_6>>>

Modal 尊重顺序并分配最优先的可用 GPU。

## 多 GPU 训练

Modal 支持单节点上的多 GPU 训练。多节点训练处于内测阶段。

### PyTorch 示例

对于重新执行入口点的框架，使用子流程或特定策略：

```python
@app.function(gpu="A100:2")
def train():
    import subprocess
    import sys
    subprocess.run(
        ["python", "train.py"],
        stdout=sys.stdout,
        stderr=sys.stderr,
        check=True,
    )
```

对于 PyTorch Lightning，将策略设置为 `ddp_spawn` 或 `ddp_notebook`。

## 性能考虑因素

**内存限制与计算限制**：
- 运行小批量模型会受到内存限制
- 较新的 GPU 的运算速度比内存访问速度更快
- 新硬件的加速可能无法证明内存受限工作负载的成本合理

**优化**：
- 尽可能使用批处理
- 在跳到 H100/B200 之前考虑 L40S
- 配置文件以识别瓶颈