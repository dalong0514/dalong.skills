<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：esm
描述：蛋白质语言模型的综合工具包，包括 ESM3（跨序列、结构和功能的生成多模式蛋白质设计）和 ESM C（高效蛋白质嵌入和表示）。在处理蛋白质序列、结构或功能预测时使用此技能；设计新型蛋白质；生成蛋白质嵌入；进行反向折叠；或进行蛋白质工程任务。支持本地模型使用和基于云的 Forge API 以进行可扩展推理。
---

# ESM：进化尺度建模

## 概述

ESM 提供最先进的蛋白质语言模型，用于理解、生成和设计蛋白质。该技能支持使用两个模型系列：用于跨序列、结构和功能的生成蛋白质设计的 ESM3，以及用于高效蛋白质表示学习和嵌入的 ESM C。

## 核心能力

### 1. 使用 ESM3 生成蛋白质序列

使用多模式生成模型生成具有所需特性的新型蛋白质序列。

**何时使用：**
- 设计具有特定功能特性的蛋白质
- 完成部分蛋白质序列
- 生成现有蛋白质的变体
- 创造具有所需结构特征的蛋白质

**基本用法：**

```python
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig

# Load model locally
model: ESM3InferenceClient = ESM3.from_pretrained("esm3-sm-open-v1").to("cuda")

# Create protein prompt
protein = ESMProtein(sequence="MPRT___KEND")  # '_' represents masked positions

# Generate completion
protein = model.generate(protein, GenerationConfig(track="sequence", num_steps=8))
print(protein.sequence)
```

**通过 Forge API 进行远程/云使用：**

<<<代码块_1>>>

请参阅 `references/esm3-api.md` 了解详细的 ESM3 模型规范、高级生成配置和多模式提示示例。

### 2.结构预测和逆折叠

使用ESM3的结构轨迹从序列或逆折叠（从结构进行序列设计）进行结构预测。

**结构预测：**

<<<代码块_2>>>

**反向折叠（结构顺序）：**

<<<代码块_3>>>

### 3. 使用 ESM C 进行蛋白质嵌入

为下游任务（例如函数预测、分类或相似性分析）生成高质量的嵌入。

**何时使用：**
- 提取蛋白质表示用于机器学习
- 计算序列相似度
- 蛋白质分类的特征提取
- 蛋白质相关任务的迁移学习

**基本用法：**

<<<代码块_4>>>

**批量处理：**

<<<代码块_5>>>

请参阅 `references/esm-c-api.md` 了解 ESM C 模型详细信息、效率比较和高级嵌入策略。

### 4. 函数条件和注释

使用 ESM3 的功能跟踪生成具有特定功能注释的蛋白质或根据序列预测功能。

**功能条件生成：**

<<<代码块_6>>>

### 5. 思想链的生成

使用 ESM3 的思想链生成方法迭代完善蛋白质设计。

```python
from esm.sdk.api import GenerationConfig

# Multi-step refinement
protein = ESMProtein(sequence="MPRT" + "_" * 100 + "KEND")

# Step 1: Generate initial structure
config = GenerationConfig(track="structure", num_steps=50)
protein = model.generate(protein, config)

# Step 2: Refine sequence based on structure
config = GenerationConfig(track="sequence", num_steps=50, temperature=0.5)
protein = model.generate(protein, config)

# Step 3: Predict function
config = GenerationConfig(track="function", num_steps=20)
protein = model.generate(protein, config)
```

### 6. 使用 Forge API 进行批处理

使用 Forge 的异步执行器有效处理多种蛋白质。

```python
from esm.sdk.forge import ESM3ForgeInferenceClient
import asyncio

client = ESM3ForgeInferenceClient(model="esm3-medium-2024-08", token="<token>")

# Async batch processing
async def batch_generate(proteins_list):
    tasks = [
        client.async_generate(protein, GenerationConfig(track="sequence"))
        for protein in proteins_list
    ]
    return await asyncio.gather(*tasks)

# Execute
proteins = [ESMProtein(sequence=f"MPRT{'_' * 50}KEND") for _ in range(10)]
results = asyncio.run(batch_generate(proteins))
```

请参阅 `references/forge-api.md` 了解详细的 Forge API 文档、身份验证、速率限制和批处理模式。

## 选型指南

**ESM3 模型（生成）：**
- `esm3-sm-open-v1` (1.4B) - 开放权重，本地使用，适合实验
- `esm3-medium-2024-08` (7B) - 质量和速度的最佳平衡（仅限 Forge）
- `esm3-large-2024-03` (98B) - 最高质量，速度较慢（仅限 Forge）

**ESM C 模型（嵌入）：**
- `esmc-300m`（30 层）- 轻量级、快速推理
- `esmc-600m`（36 层）- 平衡性能
- `esmc-6b`（80 层）- 最高表示质量

**选择标准：**
- **本地开发/测试：** 使用 `esm3-sm-open-v1` 或 `esmc-300m`
- **生产质量：** 通过 Forge 使用 `esm3-medium-2024-08`
- **最大精度：** 使用 `esm3-large-2024-03` 或 `esmc-6b`
- **高吞吐量：** 使用带有批处理执行器的 Forge API
- **成本优化：** 使用较小的模型，实施缓存策略

## 安装

**基本安装：**

```bash
uv pip install esm
```

**使用 Flash Attention（推荐用于更快的推理）：**

```bash
uv pip install esm
uv pip install flash-attn --no-build-isolation
```

**对于 Forge API 访问：**

```bash
uv pip install esm  # SDK includes Forge client
```

不需要额外的依赖项。从 https://forge.evolutionaryscale.ai 获取 Forge API 令牌

## 常见工作流程

有关详细示例和完整工作流程，请参阅 `references/workflows.md`，其中包括：
- 新颖的 GFP 设计思路
- 蛋白质变体生成和筛选
- 基于结构的序列优化
- 函数预测管道
- 基于嵌入的聚类和分析

## 参考文献

该技能包括全面的参考文档：
- `references/esm3-api.md` - ESM3模型架构、API参考、生成参数和多模态提示
- `references/esm-c-api.md` - ESM C 模型详细信息、嵌入策略和性能优化
- `references/forge-api.md` - Forge 平台文档、身份验证、批处理和部署
- `references/workflows.md` - 完整示例和常见工作流程模式

这些参考资料包含详细的 API 规范、参数说明和高级使用模式。根据特定任务的需要加载它们。

## 最佳实践

**对于生成任务：**
- 从较小的模型开始进行原型设计 (`esm3-sm-open-v1`)
- 使用温度参数来控制多样性（0.0 = 确定性，1.0 = 多样性）
- 通过复杂设计的思想链实现迭代细化
- 通过结构预测或湿实验室实验验证生成的序列

**对于嵌入任务：**
- 尽可能提高效率的批处理序列
- 用于重复分析的缓存嵌入
- 计算相似度时标准化嵌入
- 根据下游任务要求使用适当的模型大小

**对于生产部署：**
- 使用 Forge API 实现可扩展性和最新模型
- 实现 API 调用的错误处理和重试逻辑
- 监控代币使用情况并实施速率限制
- 考虑为专用基础设施部署 AWS SageMaker

## 资源和文档

- **GitHub 存储库：** https://github.com/evolutionaryscale/esm
- **锻造平台：** https://forge.evolutionaryscale.ai
- **科学论文：** Hayes 等人，《科学》(2025) - https://www.science.org/doi/10.1126/science.ads0018
- **博客文章：**
  - ESM3 发布：https://www.evolutionaryscale.ai/blog/esm3-release
  - ESM C 启动：https://www.evolutionaryscale.ai/blog/esm-cambrian
- **社区：** Slack 社区，位于 https://bit.ly/3FKwcWd
- **模型权重：** HuggingFace EvolutionaryScale 组织

## 负责任的使用

ESM 专为蛋白质工程、药物发现和科学研究中的有益应用而设计。设计新型蛋白质时遵循负责任的生物设计框架 (https://responsiblebiodesign.ai/)。在实验验证之前考虑蛋白质设计的生物安全性和伦理影响。