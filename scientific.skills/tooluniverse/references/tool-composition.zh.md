<!-- 此文件由机器翻译自 tool-composition.md -->

# ToolUniverse 中的工具组成和工作流程

## 概述

ToolUniverse 可以将多个工具链接在一起以创建复杂的科学工作流程。可以顺序或并行地组合工具来解决多步骤的研究问题。

## 顺序工具组合

按顺序执行工具，每个工具的输出都会输入下一个工具。

### 基本模式
```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

# Step 1: Get disease-associated targets
targets = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {"efoId": "EFO_0000537"}  # Hypertension
})

# Step 2: For each target, get protein structure
structures = []
for target in targets[:5]:  # First 5 targets
    structure = tu.run({
        "name": "AlphaFold_get_structure",
        "arguments": {"uniprot_id": target['uniprot_id']}
    })
    structures.append(structure)

# Step 3: Analyze structures
for structure in structures:
    analysis = tu.run({
        "name": "ProteinAnalysis_calculate_properties",
        "arguments": {"structure": structure}
    })
```

## 复杂工作流程示例

### 药物发现工作流程

从疾病到候选药物的完整工作流程：

<<<代码块_1>>>

### 基因组分析工作流程

<<<代码块_2>>>

### 临床基因组学工作流程

<<<代码块_3>>>

## 并行工具执行

当多个工具互不依赖时同时执行它们：

<<<代码块_4>>>

## 输出处理挂钩

ToolUniverse 支持后处理挂钩：
- 总结
- 文件保存
- 数据转换
- 可视化

<<<代码块_5>>>

## 最佳实践

1. **错误处理**：为工作流程中的每个工具实现 try- except 块
2. **数据验证**：在传递到下一个工具之前验证每个步骤的输出
3. **检查点**：保存长工作流程的中间结果
4. **日志**：通过复杂的工作流程跟踪进度
5. **资源管理**：考虑速率限制和计算资源
6. **模块化**：将复杂的工作流程分解为可重用的功能
7. **测试**：在组成完整的工作流程之前单独测试每个步骤