<!-- 此文件由机器翻译自 tool-discovery.md -->

# ToolUniverse 中的工具发现

## 概述

ToolUniverse 提供了多种方法来使用自然语言、关键字或嵌入来发现和搜索 600 多种科学工具。

## 发现方法

### 1. Tool_Finder（基于嵌入的搜索）

使用语义嵌入来查找相关工具。 **需要 GPU** 才能获得最佳性能。

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

# Search by natural language description
tools = tu.run({
    "name": "Tool_Finder",
    "arguments": {
        "description": "protein structure prediction",
        "limit": 10
    }
})

print(tools)
```

**何时使用：**
- 自然语言查询
- 语义相似度搜索
- 当 GPU 可用时

### 2. Tool_Finder_LLM（基于LLM的搜索）

使用 LLM 推理的基于嵌入的搜索的替代方案。 **无需 GPU**。

<<<代码块_1>>>

**何时使用：**
- 当 GPU 不可用时
- 需要推理的复杂查询
- 需要语义理解

### 3. Tool_Finder_Keyword（关键字搜索）

通过工具名称和描述进行基于关键字的快速搜索。

<<<代码块_2>>>

**何时使用：**
- 快速搜索
- 已知关键词
- 精确的术语匹配

## 列出可用工具

### 列出所有工具
<<<代码块_3>>>

### 列出有限制的工具
<<<代码块_4>>>

## 工具信息

### 获取工具详细信息
<<<代码块_5>>>

## 搜索策略

### 按域
使用特定于域的关键字：
- 生物信息学：“序列比对”、“基因组学”、“RNA-seq”
- 化学信息学：“分子动力学”、“药物设计”、“SMILES”
- 机器学习：“分类”、“预测”、“神经网络”
- 结构生物学：“蛋白质结构”、“PDB”、“晶体学”

### 按功能
按您想要实现的目标进行搜索：
- “寻找疾病基因关联”
- “预测蛋白质相互作用”
- “分析临床试验数据”
- “生成分子描述符”

### 按数据源
搜索特定数据库或 API：
- “OpenTargets”、“PubChem”、“UniProt”
-“AlphaFold”、“ChEMBL”、“PDB”
- “KEGG”、“Reactome”、“STRING”

## 最佳实践

1. **从广泛开始**：从一般术语开始，然后细化
2. **使用多种方法**：如果结果不令人满意，请尝试不同的发现方法
3. **设置适当的限制**：使用 `limit` 参数来控制结果大小（默认值：10）
4. **检查工具描述**：查看返回的工具描述以验证相关性
5. **迭代**：根据初始结果细化搜索词