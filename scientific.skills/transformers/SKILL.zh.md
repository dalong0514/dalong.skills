<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：变形金刚
描述：在使用预训练的 Transformer 模型进行自然语言处理、计算机视觉、音频或多模式任务时，应使用此技能。用于文本生成、分类、问答、翻译、摘要、图像分类、对象检测、语音识别以及自定义数据集上的模型微调。
---

# 变形金刚

## 概述

Hugging Face Transformers 库提供对数千个预训练模型的访问，这些模型适用于 NLP、计算机视觉、音频和多模态领域的任务。使用此技能来加载模型、执行推理以及对自定义数据进行微调。

## 安装

安装变压器和核心依赖项：

```bash
uv pip install torch transformers datasets evaluate accelerate
```

对于视觉任务，添加：
<<<代码块_1>>>

对于音频任务，添加：
<<<代码块_2>>>

## 身份验证

Hugging Face Hub 上的许多模型都需要身份验证。设置访问权限：

<<<代码块_3>>>

或者设置环境变量：
<<<代码块_4>>>

获取令牌：https://huggingface.co/settings/tokens

## 快速入门

使用 Pipeline API 进行快速推理，无需手动配置：

<<<代码块_5>>>

## 核心能力

### 1. 快速推理管道

用于跨许多任务进行简单、优化的推理。支持文本生成、分类、NER、问答、摘要、翻译、图像分类、对象检测、音频分类等。

**何时使用**：快速原型设计、简单的推理任务、无需自定义预处理。

有关全面的任务覆盖和优化，请参阅 `references/pipelines.md`。

### 2.模型加载与管理

加载预先训练的模型，对配置、设备放置和精度进行细粒度控制。

**何时使用**：自定义模型初始化、高级设备管理、模型检查。

有关加载模式和最佳实践，请参阅 `references/models.md`。

### 3. 文本生成

使用各种解码策略（贪婪、波束搜索、采样）和控制参数（温度、top-k、top-p）通过 LLM 生成文本。

**何时使用**：创意文本生成、代码生成、对话式 AI、文本完成。

生成策略和参数请参见`references/generation.md`。

### 4. 训练和微调

使用具有自动混合精度、分布式训练和日志记录功能的 Trainer API 对自定义数据集上的预训练模型进行微调。

**何时使用**：特定于任务的模型适应、领域适应、提高模型性能。

请参阅 `references/training.md` 了解培训工作流程和最佳实践。

### 5. 代币化

将文本转换为模型输入的标记和标记 ID，并使用填充、截断和特殊标记处理。

**何时使用**：自定义预处理管道、了解模型输入、批处理。

有关标记化的详细信息，请参阅 `references/tokenizers.md`。

## 常见模式

### 模式 1：简单推理
对于简单的任务，请使用管道：
<<<代码块_6>>>

### 模式 2：自定义模型使用
对于高级控制，分别加载模型和分词器：
```python
from transformers import AutoModelForCausalLM, AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained("model-id")
model = AutoModelForCausalLM.from_pretrained("model-id", device_map="auto")

inputs = tokenizer("text", return_tensors="pt")
outputs = model.generate(**inputs, max_new_tokens=100)
result = tokenizer.decode(outputs[0])
```

### 模式 3：微调
对于任务适应，请使用 Trainer：
```python
from transformers import Trainer, TrainingArguments

training_args = TrainingArguments(
    output_dir="./results",
    num_train_epochs=3,
    per_device_train_batch_size=8,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
)

trainer.train()
```

## 参考文档

有关特定组件的详细信息：
- **管道**：`references/pipelines.md` - 所有支持的任务和优化
- **模型**：`references/models.md` - 加载、保存和配置
- **生成**：`references/generation.md` - 文本生成策略和参数
- **训练**：`references/training.md` - 使用 Trainer API 进行微调
- **分词器**：`references/tokenizers.md` - 分词和预处理