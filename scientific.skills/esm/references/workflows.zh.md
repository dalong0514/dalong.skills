<!-- 此文件由机器翻译自 workflows.md -->

# ESM 工作流程和示例

## 概述

本文档提供了使用 ESM3 和 ESM C 的常见工作流程的完整端到端示例。每个工作流程包括设置、执行和分析代码。

## 工作流程 1：具有思路的新颖 GFP 设计

使用 ESM3 的多模式生成功能设计新型荧光蛋白。

### 目标

使用跨序列、结构和功能的思维链推理生成具有特定属性的绿色荧光蛋白 (GFP)。

### 完成实施

```python
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, GenerationConfig, FunctionAnnotation
import matplotlib.pyplot as plt

# Setup
model = ESM3.from_pretrained("esm3-sm-open-v1").to("cuda")

# Step 1: Define target properties
print("Step 1: Defining target GFP properties...")

# Create protein with desired function
target_length = 238  # Typical GFP length
protein = ESMProtein(
    sequence="_" * target_length,
    function_annotations=[
        FunctionAnnotation(
            label="green_fluorescent_protein",
            start=65,
            end=75  # Chromophore region
        )
    ]
)

# Step 2: Generate initial sequence with function conditioning
print("Step 2: Generating initial sequence...")

config = GenerationConfig(
    track="sequence",
    num_steps=target_length // 3,  # Gradual generation
    temperature=0.7  # Moderate diversity
)
protein = model.generate(protein, config)
print(f"Generated sequence: {protein.sequence[:50]}...")

# Step 3: Predict structure
print("Step 3: Predicting structure...")

config = GenerationConfig(
    track="structure",
    num_steps=target_length // 2
)
protein = model.generate(protein, config)
print(f"Structure predicted, coordinates shape: {protein.coordinates.shape}")

# Step 4: Refine sequence based on structure
print("Step 4: Refining sequence based on structure...")

# Mask regions for refinement (e.g., surface residues)
sequence_list = list(protein.sequence)
# Keep chromophore region, refine others
for i in range(0, 65):
    if i % 3 == 0:  # Refine every third position
        sequence_list[i] = '_'
for i in range(75, target_length):
    if i % 3 == 0:
        sequence_list[i] = '_'

protein.sequence = ''.join(sequence_list)

config = GenerationConfig(
    track="sequence",
    num_steps=50,
    temperature=0.5  # Lower temperature for refinement
)
protein = model.generate(protein, config)

# Step 5: Final validation
print("Step 5: Final validation...")

# Predict final structure
config = GenerationConfig(track="structure", num_steps=30)
protein = model.generate(protein, config)

# Save results
with open("novel_gfp.pdb", "w") as f:
    f.write(protein.to_pdb())

with open("novel_gfp_sequence.txt", "w") as f:
    f.write(f">Novel_GFP\n{protein.sequence}\n")

print(f"\nFinal GFP sequence:\n{protein.sequence}")
print(f"\nFunction annotations: {protein.function_annotations}")
print(f"Structure saved to: novel_gfp.pdb")
```

### 验证步骤

<<<代码块_1>>>

## 工作流程 2：蛋白质变体库生成

生成并分析用于定向进化的蛋白质变体库。

### 目标

通过定向诱变创建亲本蛋白的变体，同时保持结构完整性。

### 完成实施

<<<代码块_2>>>

## 工作流程 3：基于结构的序列优化

优化蛋白质序列以提高稳定性，同时保持功能。

### 目标

给定蛋白质结构，设计保持折叠但具有改进特性的序列。

### 完成实施

<<<代码块_3>>>

## 工作流程 4：函数预测管道

使用 ESM3 和 ESM C 根据序列预测蛋白质功能。

### 目标

使用生成 (ESM3) 和嵌入 (ESM C) 方法构建预测蛋白质功能的管道。

### 完成实施

<<<代码块_4>>>

## 工作流程 5：基于嵌入的聚类和分析

使用 ESM C 嵌入对大型蛋白质数据集进行聚类和分析。

### 完成实施

<<<代码块_5>>>

## 其他工作流程提示

### 大型数据集的内存管理

<<<代码块_6>>>

### 并行处理

```python
from concurrent.futures import ThreadPoolExecutor
import asyncio

def parallel_workflow(sequences, n_workers=4):
    """Process sequences in parallel."""

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(process_sequence, sequences))

    return results
```

这些工作流程为常见 ESM 用例提供了全面的示例。使它们适应您的特定需求，并始终通过适当的生物实验验证结果。