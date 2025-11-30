<!-- 此文件由机器翻译自 data_module.md -->

# LightningDataModule - 综合指南

## 概述

LightningDataModule 是一个可重用、可共享的类，封装了 PyTorch Lightning 中的所有数据处理步骤。它通过标准化数据集的管理和跨项目共享的方式解决了分散的数据准备逻辑问题。

## 它解决的核心问题

在传统的 PyTorch 工作流程中，数据处理分散在多个文件中，因此很难回答以下问题：
- “你用了什么分裂？”
- “应用了哪些变换？”
- “数据是如何准备的？”

DataModules 集中这些信息以实现可重复性和可重用性。

## 五个处理步骤

DataModule 将数据处理分为五个阶段：

1. **下载/标记化/处理** - 初始数据采集
2. **清理并保存** - 将处理后的数据保存到磁盘
3. **加载到数据集** - 创建 PyTorch 数据集对象
4. **应用转换** - 数据增强、标准化等。
5. **Wrap in DataLoader** - 配置批处理和加载

## 主要方法

### `prepare_data()`
下载并处理数据。仅在单个进程上运行一次（非分布式）。

**用于：**
- 下载数据集
- 标记文本
- 将处理后的数据保存到磁盘

**重要提示：** 不要在此处设置状态（例如 self.x = y）。状态不会转移到其他进程。

**示例：**
```python
def prepare_data(self):
    # Download data (runs once)
    download_dataset("http://example.com/data.zip", "data/")

    # Tokenize and save (runs once)
    tokenize_and_save("data/raw/", "data/processed/")
```

### `setup(stage)`
创建数据集并应用转换。在分布式训练的每个进程上运行。

**参数：**
- `stage` - '适合'、'验证'、'测试'或'预测'

**用于：**
- 创建训练/验证/测试分割
- 构建数据集对象
- 应用变换
- 设置状态（self.train_dataset = ...）

**示例：**
<<<代码块_1>>>

### `train_dataloader()`
返回训练数据加载器。

**示例：**
<<<代码块_2>>>

### `val_dataloader()`
返回验证数据加载器。

**示例：**
<<<代码块_3>>>

### `test_dataloader()`
返回测试数据加载器。

**示例：**
<<<代码块_4>>>

### `predict_dataloader()`
返回预测数据加载器。

**示例：**
<<<代码块_5>>>

## 完整示例

<<<代码块_6>>>

## 用法

```python
# Create DataModule
dm = MyDataModule(data_dir="./data", batch_size=64, num_workers=8)

# Use with Trainer
trainer = L.Trainer(max_epochs=10)
trainer.fit(model, datamodule=dm)

# Test
trainer.test(model, datamodule=dm)

# Predict
predictions = trainer.predict(model, datamodule=dm)

# Or use standalone in PyTorch
dm.prepare_data()
dm.setup(stage='fit')
train_loader = dm.train_dataloader()

for batch in train_loader:
    # Your training code
    pass
```

## 附加挂钩

### `transfer_batch_to_device(batch, device, dataloader_idx)`
用于将批次移动到设备的自定义逻辑。

**示例：**
```python
def transfer_batch_to_device(self, batch, device, dataloader_idx):
    # Custom transfer logic
    if isinstance(batch, dict):
        return {k: v.to(device) for k, v in batch.items()}
    return super().transfer_batch_to_device(batch, device, dataloader_idx)
```

### `on_before_batch_transfer(batch, dataloader_idx)`
在传输到设备之前增加或修改批次（在 CPU 上运行）。

**示例：**
```python
def on_before_batch_transfer(self, batch, dataloader_idx):
    # Apply CPU-based augmentations
    batch['image'] = apply_augmentation(batch['image'])
    return batch
```

### `on_after_batch_transfer(batch, dataloader_idx)`
传输到设备后增加或修改批次（在 GPU 上运行）。

**示例：**
```python
def on_after_batch_transfer(self, batch, dataloader_idx):
    # Apply GPU-based augmentations
    batch['image'] = gpu_augmentation(batch['image'])
    return batch
```

### `state_dict()` / `load_state_dict(state_dict)`
保存和恢复 DataModule 状态以进行检查点。

**示例：**
```python
def state_dict(self):
    return {"current_fold": self.current_fold}

def load_state_dict(self, state_dict):
    self.current_fold = state_dict["current_fold"]
```

### `teardown(stage)`
训练/测试/预测后的清理操作。

**示例：**
```python
def teardown(self, stage):
    # Clean up resources
    if stage == 'fit':
        self.train_dataset = None
        self.val_dataset = None
```

## 高级模式

### 多个验证/测试数据加载器

返回 DataLoaders 的列表或字典：

```python
def val_dataloader(self):
    return [
        DataLoader(self.val_dataset_1, batch_size=32),
        DataLoader(self.val_dataset_2, batch_size=32)
    ]

# Or with names (for logging)
def val_dataloader(self):
    return {
        "val_easy": DataLoader(self.val_easy, batch_size=32),
        "val_hard": DataLoader(self.val_hard, batch_size=32)
    }

# In LightningModule
def validation_step(self, batch, batch_idx, dataloader_idx=0):
    if dataloader_idx == 0:
        # Handle val_dataset_1
        pass
    else:
        # Handle val_dataset_2
        pass
```

### 交叉验证

```python
class CrossValidationDataModule(L.LightningDataModule):
    def __init__(self, data_dir, batch_size, num_folds=5):
        super().__init__()
        self.data_dir = data_dir
        self.batch_size = batch_size
        self.num_folds = num_folds
        self.current_fold = 0

    def setup(self, stage=None):
        full_dataset = MyDataset(self.data_dir)
        fold_size = len(full_dataset) // self.num_folds

        # Create fold indices
        indices = list(range(len(full_dataset)))
        val_start = self.current_fold * fold_size
        val_end = val_start + fold_size

        val_indices = indices[val_start:val_end]
        train_indices = indices[:val_start] + indices[val_end:]

        self.train_dataset = Subset(full_dataset, train_indices)
        self.val_dataset = Subset(full_dataset, val_indices)

    def set_fold(self, fold):
        self.current_fold = fold

    def state_dict(self):
        return {"current_fold": self.current_fold}

    def load_state_dict(self, state_dict):
        self.current_fold = state_dict["current_fold"]

# Usage
dm = CrossValidationDataModule("./data", batch_size=32, num_folds=5)

for fold in range(5):
    dm.set_fold(fold)
    trainer = L.Trainer(max_epochs=10)
    trainer.fit(model, datamodule=dm)
```

### 超参数保存

```python
class MyDataModule(L.LightningDataModule):
    def __init__(self, data_dir, batch_size=32, num_workers=4):
        super().__init__()
        # Save hyperparameters
        self.save_hyperparameters()

    def setup(self, stage=None):
        # Access via self.hparams
        print(f"Batch size: {self.hparams.batch_size}")
```

## 最佳实践

### 1. 分离prepare_data和setup
- `prepare_data()` - 下载/进程（单个进程，无状态）
- `setup()` - 创建数据集（每个进程，设置状态）

### 2.使用stage参数
检查 `setup()` 中的阶段以避免不必要的工作：

```python
def setup(self, stage):
    if stage == 'fit':
        # Only load train/val data when fitting
        self.train_dataset = ...
        self.val_dataset = ...
    elif stage == 'test':
        # Only load test data when testing
        self.test_dataset = ...
```

### 3. 用于 GPU 训练的 Pin 内存
在 DataLoaders 中启用 `pin_memory=True` 以加快 GPU 传输速度：

```python
def train_dataloader(self):
    return DataLoader(..., pin_memory=True)
```

### 4. 使用持久的员工
防止工作进程在纪元之间重新启动：

```python
def train_dataloader(self):
    return DataLoader(
        ...,
        num_workers=4,
        persistent_workers=True
    )
```

### 5. 避免验证/测试中的随机播放
切勿打乱验证或测试数据：

```python
def val_dataloader(self):
    return DataLoader(..., shuffle=False)  # Never True
```

### 6. 使 DataModule 可重用
接受`__init__`中的配置参数：

```python
class MyDataModule(L.LightningDataModule):
    def __init__(self, data_dir, batch_size, num_workers, augment=True):
        super().__init__()
        self.save_hyperparameters()
```

### 7. 文档数据结构
添加解释数据格式和期望的文档字符串：

```python
class MyDataModule(L.LightningDataModule):
    """
    DataModule for XYZ dataset.

    Data format: (image, label) tuples
    - image: torch.Tensor of shape (C, H, W)
    - label: int in range [0, num_classes)

    Args:
        data_dir: Path to data directory
        batch_size: Batch size for dataloaders
        num_workers: Number of data loading workers
    """
```

## 常见陷阱

### 1.在prepare_data中设置状态
**错误：**
```python
def prepare_data(self):
    self.dataset = load_data()  # State not transferred to other processes!
```

**正确：**
```python
def prepare_data(self):
    download_data()  # Only download, no state

def setup(self, stage):
    self.dataset = load_data()  # Set state here
```

### 2. 不使用 stage 参数
**效率低下：**
```python
def setup(self, stage):
    self.train_dataset = load_train()
    self.val_dataset = load_val()
    self.test_dataset = load_test()  # Loads even when just fitting
```

**高效：**
```python
def setup(self, stage):
    if stage == 'fit':
        self.train_dataset = load_train()
        self.val_dataset = load_val()
    elif stage == 'test':
        self.test_dataset = load_test()
```

### 3. 忘记返回 DataLoaders
**错误：**
```python
def train_dataloader(self):
    DataLoader(self.train_dataset, ...)  # Forgot return!
```

**正确：**
```python
def train_dataloader(self):
    return DataLoader(self.train_dataset, ...)
```

## 与训练器集成

```python
# Initialize DataModule
dm = MyDataModule(data_dir="./data", batch_size=64)

# All data loading is handled by DataModule
trainer = L.Trainer(max_epochs=10)
trainer.fit(model, datamodule=dm)

# DataModule handles validation too
trainer.validate(model, datamodule=dm)

# And testing
trainer.test(model, datamodule=dm)

# And prediction
predictions = trainer.predict(model, datamodule=dm)
```