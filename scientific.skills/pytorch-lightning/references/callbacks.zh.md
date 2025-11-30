<!-- 此文件由机器翻译自 callbacks.md -->

# 回调 - 综合指南

## 概述

回调可以在训练中添加任意独立的程序，而不会扰乱您的 LightningModule 研究代码。他们在训练生命周期中在特定的钩子上执行自定义逻辑。

## 架构

Lightning 通过三个组件组织训练逻辑：
- **培训师** - 工程基础设施
- **LightningModule** - 研究代码
- **回调** - 非必要功能（监控、检查点、自定义行为）

## 创建自定义回调

基本结构：

```python
from lightning.pytorch.callbacks import Callback

class MyCustomCallback(Callback):
    def on_train_start(self, trainer, pl_module):
        print("Training is starting!")

    def on_train_end(self, trainer, pl_module):
        print("Training is done!")

# Use with Trainer
trainer = L.Trainer(callbacks=[MyCustomCallback()])
```

## 内置回调

### 模型检查点

根据监控的指标保存模型。

**关键参数：**
- `dirpath` - 保存检查点的目录
- `filename` - 检查点文件名模式
- `monitor` - 要监控的指标
- `mode` - 受监控指标的“min”或“max”
- `save_top_k` - 要保留的最佳模型数量
- `save_last` - 保存最后一个纪元检查点
- `every_n_epochs` - 每 N 个时期保存一次
- `save_on_train_epoch_end` - 在训练纪元结束与验证结束时保存

**示例：**
<<<代码块_1>>>

**访问保存的检查点：**
<<<代码块_2>>>

### 提前停止

当监控指标停止改善时停止训练。

**关键参数：**
- `monitor` - 要监控的指标
- `patience` - 训练停止后没有改善的时期数
- `mode` - 受监控指标的“min”或“max”
- `min_delta` - 符合改进资格的最小更改
- `verbose` - 打印消息
- `strict` - 如果未找到受监控的指标，则会崩溃

**示例：**
<<<代码块_3>>>

### 学习率监视器

跟踪调度程序的学习率变化。

**关键参数：**
- `logging_interval` - 何时记录：“step”或“epoch”
- `log_momentum` - 还记录动量值

**示例：**
<<<代码块_4>>>

### 设备统计监视器

记录设备性能指标（GPU/CPU/TPU）。

**关键参数：**
- `cpu_stats` - 记录 CPU 统计信息

**示例：**
<<<代码块_5>>>

### ModelSummary / RichModelSummary

显示模型架构和参数计数。

**示例：**
<<<代码块_6>>>

### 计时器

跟踪并限制培训持续时间。

**关键参数：**
- `duration` - 最大训练时间（timedelta 或 dict）
- `interval` - 检查间隔：“step”、“epoch”或“batch”

**示例：**
```python
from lightning.pytorch.callbacks import Timer
from datetime import timedelta

# Limit training to 1 hour
timer = Timer(duration=timedelta(hours=1))

# Or using dict
timer = Timer(duration={"hours": 23, "minutes": 30})

trainer = L.Trainer(callbacks=[timer])
```

### 批量大小查找器

自动找到最佳批量大小。

**示例：**
```python
from lightning.pytorch.callbacks import BatchSizeFinder

batch_finder = BatchSizeFinder(mode="power", steps_per_trial=3)

trainer = L.Trainer(callbacks=[batch_finder])
trainer.fit(model, datamodule=dm)

# Optimal batch size is set automatically
```

### 梯度累积调度器

动态安排梯度累积步骤。

**示例：**
```python
from lightning.pytorch.callbacks import GradientAccumulationScheduler

# Accumulate 4 batches for first 5 epochs, then 2 batches
accumulator = GradientAccumulationScheduler(scheduling={0: 4, 5: 2})

trainer = L.Trainer(callbacks=[accumulator])
```

### 随机权重平均 (SWA)

应用随机权重平均以获得更好的泛化能力。

**示例：**
```python
from lightning.pytorch.callbacks import StochasticWeightAveraging

swa = StochasticWeightAveraging(swa_lrs=1e-2, swa_epoch_start=0.8)

trainer = L.Trainer(callbacks=[swa])
```

## 自定义回调示例

### 简单的日志回调

```python
class MetricsLogger(Callback):
    def __init__(self):
        self.metrics = []

    def on_validation_end(self, trainer, pl_module):
        # Access logged metrics
        metrics = trainer.callback_metrics
        self.metrics.append(dict(metrics))
        print(f"Validation metrics: {metrics}")
```

### 梯度监控回调

```python
class GradientMonitor(Callback):
    def on_after_backward(self, trainer, pl_module):
        # Log gradient norms
        for name, param in pl_module.named_parameters():
            if param.grad is not None:
                grad_norm = param.grad.norm().item()
                pl_module.log(f"grad_norm/{name}", grad_norm)
```

### 自定义检查点回调

```python
class CustomCheckpoint(Callback):
    def __init__(self, save_dir):
        self.save_dir = save_dir

    def on_train_epoch_end(self, trainer, pl_module):
        epoch = trainer.current_epoch
        if epoch % 5 == 0:  # Save every 5 epochs
            filepath = f"{self.save_dir}/custom-{epoch}.ckpt"
            trainer.save_checkpoint(filepath)
            print(f"Saved checkpoint: {filepath}")
```

### 模型冻结回调

```python
class FreezeUnfreeze(Callback):
    def __init__(self, freeze_until_epoch=10):
        self.freeze_until_epoch = freeze_until_epoch

    def on_train_epoch_start(self, trainer, pl_module):
        epoch = trainer.current_epoch

        if epoch < self.freeze_until_epoch:
            # Freeze backbone
            for param in pl_module.backbone.parameters():
                param.requires_grad = False
        else:
            # Unfreeze backbone
            for param in pl_module.backbone.parameters():
                param.requires_grad = True
```

### 学习率查找器回调

```python
class LRFinder(Callback):
    def __init__(self, min_lr=1e-5, max_lr=1e-1, num_steps=100):
        self.min_lr = min_lr
        self.max_lr = max_lr
        self.num_steps = num_steps
        self.lrs = []
        self.losses = []

    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        if batch_idx >= self.num_steps:
            trainer.should_stop = True
            return

        # Exponential LR schedule
        lr = self.min_lr * (self.max_lr / self.min_lr) ** (batch_idx / self.num_steps)
        optimizer = trainer.optimizers[0]
        for param_group in optimizer.param_groups:
            param_group['lr'] = lr

        self.lrs.append(lr)
        self.losses.append(outputs['loss'].item())

    def on_train_end(self, trainer, pl_module):
        # Plot LR vs Loss
        import matplotlib.pyplot as plt
        plt.plot(self.lrs, self.losses)
        plt.xscale('log')
        plt.xlabel('Learning Rate')
        plt.ylabel('Loss')
        plt.savefig('lr_finder.png')
```

### 预测保存器回调

```python
class PredictionSaver(Callback):
    def __init__(self, save_path):
        self.save_path = save_path
        self.predictions = []

    def on_predict_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        self.predictions.append(outputs)

    def on_predict_end(self, trainer, pl_module):
        # Save all predictions
        torch.save(self.predictions, self.save_path)
        print(f"Predictions saved to {self.save_path}")
```

## 可用的钩子

### 设置和拆卸
- `setup(trainer, pl_module, stage)` - 在拟合/测试/预测开始时调用
- `teardown(trainer, pl_module, stage)` - 在拟合/测试/预测结束时调用

### 训练生命周期
- `on_fit_start(trainer, pl_module)` - 在拟合开始时调用
- `on_fit_end(trainer, pl_module)` - 在拟合结束时调用
- `on_train_start(trainer, pl_module)` - 在训练开始时调用
- `on_train_end(trainer, pl_module)` - 在训练结束时调用

### 纪元边界
- `on_train_epoch_start(trainer, pl_module)` - 在训练纪元开始时调用
- `on_train_epoch_end(trainer, pl_module)` - 在训练时期结束时调用
- `on_validation_epoch_start(trainer, pl_module)` - 在验证开始时调用
- `on_validation_epoch_end(trainer, pl_module)` - 在验证结束时调用
- `on_test_epoch_start(trainer, pl_module)` - 在测试开始时调用
- `on_test_epoch_end(trainer, pl_module)` - 在测试结束时调用

### 批次边界
- `on_train_batch_start(trainer, pl_module, batch, batch_idx)` - 训练批次之前
- `on_train_batch_end(trainer, pl_module, outputs, batch, batch_idx)` - 训练批次后
- `on_validation_batch_start(trainer, pl_module, batch, batch_idx)` - 验证批次之前
- `on_validation_batch_end(trainer, pl_module, outputs, batch, batch_idx)` - 验证批次后

### 渐变事件
- `on_before_backward(trainer, pl_module, loss)` - 在loss.backward()之前
- `on_after_backward(trainer, pl_module)` - 在loss.backward()之后
- `on_before_optimizer_step(trainer, pl_module, optimizer)` - 在optimizer.step()之前

### 检查站活动
- `on_save_checkpoint(trainer, pl_module, checkpoint)` - 保存检查点时
- `on_load_checkpoint(trainer, pl_module, checkpoint)` - 加载检查点时

### 异常处理
- `on_exception(trainer, pl_module, exception)` - 发生异常时

## 状态管理
对于需要跨检查点持久化的回调：

```python
class StatefulCallback(Callback):
    def __init__(self):
        self.counter = 0

    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        self.counter += 1

    def state_dict(self):
        return {"counter": self.counter}

    def load_state_dict(self, state_dict):
        self.counter = state_dict["counter"]

    @property
    def state_key(self):
        # Unique identifier for this callback
        return "my_stateful_callback"
```

## 最佳实践

### 1. 保持回调隔离
每个回调应该是自包含且独立的：

```python
# Good: Self-contained
class MyCallback(Callback):
    def __init__(self):
        self.data = []

    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        self.data.append(outputs['loss'].item())

# Bad: Depends on external state
global_data = []

class BadCallback(Callback):
    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        global_data.append(outputs['loss'].item())  # External dependency
```

### 2.避免回调间依赖
回调不应依赖于其他回调：

```python
# Bad: Callback B depends on Callback A
class CallbackA(Callback):
    def __init__(self):
        self.value = 0

class CallbackB(Callback):
    def __init__(self, callback_a):
        self.callback_a = callback_a  # Tight coupling

# Good: Independent callbacks
class CallbackA(Callback):
    def __init__(self):
        self.value = 0

class CallbackB(Callback):
    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        # Access trainer state instead
        value = trainer.callback_metrics.get('metric')
```

### 3.永远不要手动调用回调方法
让 Lightning 自动调用回调：

```python
# Bad: Manual invocation
callback = MyCallback()
callback.on_train_start(trainer, model)  # Don't do this

# Good: Let Trainer handle it
trainer = L.Trainer(callbacks=[MyCallback()])
```

### 4. 任何执行订单的设计
回调可以按任何顺序执行，因此不要依赖特定的顺序：

```python
# Good: Order-independent
class GoodCallback(Callback):
    def on_train_epoch_end(self, trainer, pl_module):
        # Use trainer state, not other callbacks
        metrics = trainer.callback_metrics
        self.log_metrics(metrics)
```

### 5. 对非必要逻辑使用回调
将核心研究代码保留在LightningModule中，使用回调来实现辅助功能：

```python
# Good separation
class MyModel(L.LightningModule):
    # Core research logic here
    def training_step(self, batch, batch_idx):
        return loss

# Non-essential monitoring in callback
class MonitorCallback(Callback):
    def on_validation_end(self, trainer, pl_module):
        # Monitoring logic
        pass
```

## 常见模式

### 组合多个回调

```python
from lightning.pytorch.callbacks import (
    ModelCheckpoint,
    EarlyStopping,
    LearningRateMonitor,
    DeviceStatsMonitor
)

callbacks = [
    ModelCheckpoint(monitor="val_loss", mode="min", save_top_k=3),
    EarlyStopping(monitor="val_loss", patience=10, mode="min"),
    LearningRateMonitor(logging_interval="step"),
    DeviceStatsMonitor()
]

trainer = L.Trainer(callbacks=callbacks)
```

### 有条件回调激活

```python
class ConditionalCallback(Callback):
    def __init__(self, activate_after_epoch=10):
        self.activate_after_epoch = activate_after_epoch

    def on_train_epoch_end(self, trainer, pl_module):
        if trainer.current_epoch >= self.activate_after_epoch:
            # Only active after specified epoch
            self.do_something(trainer, pl_module)
```

### 多阶段训练回调

```python
class MultiStageTraining(Callback):
    def __init__(self, stage_epochs=[10, 20, 30]):
        self.stage_epochs = stage_epochs
        self.current_stage = 0

    def on_train_epoch_start(self, trainer, pl_module):
        epoch = trainer.current_epoch

        if epoch in self.stage_epochs:
            self.current_stage += 1
            print(f"Entering stage {self.current_stage}")

            # Adjust learning rate for new stage
            for optimizer in trainer.optimizers:
                for param_group in optimizer.param_groups:
                    param_group['lr'] *= 0.1
```