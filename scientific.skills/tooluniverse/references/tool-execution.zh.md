<!-- 此文件由机器翻译自 tool-execution.md -->

# ToolUniverse 中的工具执行

## 概述

使用 `run()` 方法通过 ToolUniverse 的标准化接口执行各个工具。

## 基本工具执行

### 标准模式
```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

# Execute a tool
result = tu.run({
    "name": "tool_name_here",
    "arguments": {
        "param1": "value1",
        "param2": "value2"
    }
})

print(result)
```

## 现实世界的例子

### 示例 1：疾病-目标关联 (OpenTargets)
<<<代码块_1>>>

### 示例 2：蛋白质结构预测
<<<代码块_2>>>

### 示例 3：化学性质计算
<<<代码块_3>>>

### 示例 4：基因表达分析
<<<代码块_4>>>

## 工具执行工作流程

### 1. 发现该工具
<<<代码块_5>>>

### 2. 检查工具参数
<<<代码块_6>>>

### 3. 使用正确的参数执行
```python
# Execute the tool
result = tu.run({
    "name": "KEGG_pathway_enrichment",
    "arguments": {
        "gene_list": ["TP53", "BRCA1", "EGFR"],
        "organism": "hsa"  # Homo sapiens
    }
})
```

## 处理工具结果

### 检查结果类型
```python
result = tu.run({
    "name": "some_tool",
    "arguments": {"param": "value"}
})

# Results can be various types
if isinstance(result, dict):
    print("Dictionary result")
elif isinstance(result, list):
    print(f"List with {len(result)} items")
elif isinstance(result, str):
    print("String result")
```

### 处理结果
```python
# Example: Processing multiple results
results = tu.run({
    "name": "PubMed_search",
    "arguments": {
        "query": "cancer immunotherapy",
        "max_results": 10
    }
})

for idx, paper in enumerate(results, 1):
    print(f"{idx}. {paper['title']}")
    print(f"   PMID: {paper['pmid']}")
    print(f"   Authors: {', '.join(paper['authors'][:3])}")
    print()
```

## 错误处理

```python
try:
    result = tu.run({
        "name": "some_tool",
        "arguments": {"param": "value"}
    })
except Exception as e:
    print(f"Tool execution failed: {e}")
    # Check if tool exists
    # Verify parameter names and types
    # Review tool documentation
```

## 最佳实践

1. **验证工具参数**：执行前始终检查所需参数
2. **从简单开始**：在复杂的工作流程之前使用简单的案例进行测试
3. **适当处理结果**：检查结果类型和结构
4. **错误恢复**：为生产代码实现 try- except 块
5. **文档**：查看参数要求和输出格式的工具描述
6. **速率限制**：了解远程工具的 API 速率限制
7. **数据验证**：验证输入数据格式（例如 SMILES、UniProt ID、基因符号）