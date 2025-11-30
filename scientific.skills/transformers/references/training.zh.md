<!-- 此文件由机器翻译自 training.md -->

# 训练和微调

## 概述

使用 Trainer API 在自定义数据集上微调预训练模型。 Trainer 处理训练循环、梯度累积、混合精度、日志记录和检查点。

## 基本微调工作流程

### 第 1 步：加载和预处理数据

```python
from datasets import load_dataset

# Load dataset
dataset = load_dataset("yelp_review_full")
train_dataset = dataset["train"]
eval_dataset = dataset["test"]

# Tokenize
from transformers import AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

def tokenize_function(examples):
    return tokenizer(
        examples["text"],
        padding="max_length",
        truncation=True,
        max_length=512
    )

train_dataset = train_dataset.map(tokenize_function, batched=True)
eval_dataset = eval_dataset.map(tokenize_function, batched=True)
```

### 第 2 步：加载模型

<<<代码块_1>>>

### 步骤 3：定义指标

<<<代码块_2>>>

### 步骤 4：配置训练

<<<代码块_3>>>

### 步骤 5：创建训练器并训练

<<<代码块_4>>>

### 第 6 步：保存模型

<<<代码块_5>>>

## TrainingArguments 参数

### 基本参数

**output_dir**：检查点和日志的目录
<<<代码块_6>>>

**num_train_epochs**：训练纪元数
```python
num_train_epochs=3
```

**per_device_train_batch_size**：每个 GPU/CPU 的批量大小
```python
per_device_train_batch_size=8
```

**learning_rate**：优化器学习率
```python
learning_rate=2e-5  # Common for BERT-style models
learning_rate=5e-5  # Common for smaller models
```

**weight_decay**：L2 正则化
```python
weight_decay=0.01
```

### 评估和保存

**eval_strategy**：何时评估（“no”、“steps”、“epoch”）
```python
eval_strategy="epoch"  # Evaluate after each epoch
eval_strategy="steps"  # Evaluate every eval_steps
```

**save_strategy**：何时保存检查点
```python
save_strategy="epoch"
save_strategy="steps"
save_steps=500
```

**load_best_model_at_end**：训练后加载最佳检查点
```python
load_best_model_at_end=True
metric_for_best_model="accuracy"  # Metric to compare
```

### 优化

**gradient_accumulation_steps**：累积多个步骤的梯度
```python
gradient_accumulation_steps=4  # Effective batch size = batch_size * 4
```

**fp16**：启用混合精度（NVIDIA GPU）
```python
fp16=True
```

**bf16**：启用 bfloat16（较新的 GPU）
```python
bf16=True
```

**gradient_checkpointing**：用计算换取内存
```python
gradient_checkpointing=True  # Slower but uses less memory
```

**optim**：优化器选择
```python
optim="adamw_torch"  # Default
optim="adamw_8bit"    # 8-bit Adam (requires bitsandbytes)
optim="adafactor"     # Memory-efficient alternative
```

### 学习率调度

**lr_scheduler_type**：学习率调度
```python
lr_scheduler_type="linear"       # Linear decay
lr_scheduler_type="cosine"       # Cosine annealing
lr_scheduler_type="constant"     # No decay
lr_scheduler_type="constant_with_warmup"
```

**warmup_steps** 或 **warmup_ratio**：热身期
```python
warmup_steps=500
# Or
warmup_ratio=0.1  # 10% of total steps
```

### 日志记录

**logging_dir**：TensorBoard 日志目录
```python
logging_dir="./logs"
```

**logging_steps**：每 N 步骤记录一次
```python
logging_steps=10
```

**report_to**：记录集成
```python
report_to=["tensorboard"]
report_to=["wandb"]
report_to=["tensorboard", "wandb"]
```

### 分布式训练

**ddp_backend**：分布式后端
```python
ddp_backend="nccl"  # For multi-GPU
```

**deepspeed**：DeepSpeed 配置文件
```python
deepspeed="ds_config.json"
```

## 数据收集器

处理动态填充和特殊预处理：

### DataCollatorWithPadding

将序列填充到批次中最长的序列：
```python
from transformers import DataCollatorWithPadding

data_collator = DataCollatorWithPadding(tokenizer=tokenizer)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    data_collator=data_collator,
)
```

### DataCollatorForLanguageModeling

对于掩码语言建模：
```python
from transformers import DataCollatorForLanguageModeling

data_collator = DataCollatorForLanguageModeling(
    tokenizer=tokenizer,
    mlm=True,
    mlm_probability=0.15
)
```

### DataCollatorForSeq2Seq

对于序列到序列任务：
```python
from transformers import DataCollatorForSeq2Seq

data_collator = DataCollatorForSeq2Seq(
    tokenizer=tokenizer,
    model=model,
    padding=True
)
```

## 定制培训

### 定制培训师

重写自定义行为的方法：

```python
from transformers import Trainer

class CustomTrainer(Trainer):
    def compute_loss(self, model, inputs, return_outputs=False):
        labels = inputs.pop("labels")
        outputs = model(**inputs)
        logits = outputs.logits

        # Custom loss computation
        loss_fct = torch.nn.CrossEntropyLoss(weight=class_weights)
        loss = loss_fct(logits.view(-1, self.model.config.num_labels), labels.view(-1))

        return (loss, outputs) if return_outputs else loss
```

### 自定义回调

监控培训：

```python
from transformers import TrainerCallback

class CustomCallback(TrainerCallback):
    def on_epoch_end(self, args, state, control, **kwargs):
        print(f"Epoch {state.epoch} completed")
        # Custom logic here
        return control

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    callbacks=[CustomCallback],
)
```

## 高级训练技巧

### 参数高效微调 (PEFT)

使用 LoRA 进行高效微调：

```python
from peft import LoraConfig, get_peft_model

lora_config = LoraConfig(
    r=16,
    lora_alpha=32,
    target_modules=["query", "value"],
    lora_dropout=0.05,
    bias="none",
    task_type="SEQ_CLS"
)

model = get_peft_model(model, lora_config)
model.print_trainable_parameters()  # Shows reduced parameter count

# Train normally with Trainer
trainer = Trainer(model=model, args=training_args, ...)
trainer.train()
```

### 梯度检查点

以速度为代价减少内存：

```python
model.gradient_checkpointing_enable()

training_args = TrainingArguments(
    gradient_checkpointing=True,
    ...
)
```

### 混合精度训练

```python
training_args = TrainingArguments(
    fp16=True,  # For NVIDIA GPUs with Tensor Cores
    # or
    bf16=True,  # For newer GPUs (A100, H100)
    ...
)
```

### DeepSpeed 集成

对于非常大的模型：

```python
# ds_config.json
{
  "train_batch_size": 16,
  "gradient_accumulation_steps": 1,
  "optimizer": {
    "type": "AdamW",
    "params": {
      "lr": 2e-5
    }
  },
  "fp16": {
    "enabled": true
  },
  "zero_optimization": {
    "stage": 2
  }
}
```

```python
training_args = TrainingArguments(
    deepspeed="ds_config.json",
    ...
)
```

## 训练技巧

### 超参数调整

共同的出发点：
- **学习率**：类似 BERT 的模型为 2e-5 至 5e-5，较小模型为 1e-4 至 1e-3
- **批量大小**：8-32，取决于 GPU 内存
- **Epochs**：2-4 用于微调，更多用于域适应
- **热身**：总步数的 10%

使用 Optuna 进行超参数搜索：

```python
def model_init():
    return AutoModelForSequenceClassification.from_pretrained(
        "bert-base-uncased",
        num_labels=5
    )

def optuna_hp_space(trial):
    return {
        "learning_rate": trial.suggest_float("learning_rate", 1e-5, 5e-5, log=True),
        "per_device_train_batch_size": trial.suggest_categorical("per_device_train_batch_size", [8, 16, 32]),
        "num_train_epochs": trial.suggest_int("num_train_epochs", 2, 5),
    }

trainer = Trainer(model_init=model_init, args=training_args, ...)
best_trial = trainer.hyperparameter_search(
    direction="maximize",
    backend="optuna",
    hp_space=optuna_hp_space,
    n_trials=10,
)
```

### 监控培训

使用张量板：
```bash
tensorboard --logdir ./logs
```

或者权重和偏差：
```python
import wandb
wandb.init(project="my-project")

training_args = TrainingArguments(
    report_to=["wandb"],
    ...
)
```

### 恢复培训

从检查点恢复：
```python
trainer.train(resume_from_checkpoint="./results/checkpoint-1000")
```

## 常见问题

**CUDA内存不足：**
- 减少批量
- 启用梯度检查点
- 使用梯度累积
- 使用8位优化器

**过度拟合：**
- 增加权重衰减
- 添加辍学
- 使用提前停止
- 减少模型大小或训练周期

**慢速训练：**
- 增加批量大小
- 启用混合精度（fp16/bf16）
- 使用多个GPU
- 优化数据加载

## 最佳实践

1. **从小处开始**：首先在小数据集子集上进行测试
2. **使用评估**：监控验证指标
3. **保存检查点**：启用save_strategy
4. **广泛记录**：使用 TensorBoard 或 W&B
5. **尝试不同的学习率**：从2e-5开始
6. **使用热身**：帮助训练稳定性
7. **启用混合精度**：更快的训练
8. **考虑 PEFT**：适用于资源有限的大型模型