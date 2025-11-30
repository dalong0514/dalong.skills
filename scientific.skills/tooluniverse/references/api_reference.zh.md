<!-- 此文件由机器翻译自 api_reference.md -->

# ToolUniverse Python API 参考

## 核心课程

### 工具宇宙

用于与 ToolUniverse 生态系统交互的主类。

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
```

#### 方法

##### `load_tools()`
将所有可用工具加载到 ToolUniverse 实例中。

<<<代码块_1>>>

**返回：** 无

**副作用：** 将 600 多个工具加载到内存中以供发现和执行。

---

##### `run(tool_config)`
使用指定参数执行工具。

**参数：**
- `tool_config` (dict)：带有键的配置字典：
  - `name` (str)：要执行的工具名称
  - `arguments` (dict)：特定于工具的参数

**返回：** 工具特定的输出（dict、list、str 或其他类型）

**示例：**
<<<代码块_2>>>

---

##### `list_tools(limit=None)`
列出所有可用工具或子集。

**参数：**
- `limit`（int，可选）：要返回的工具的最大数量。如果没有，则返回所有工具。

**返回：** 工具字典列表

**示例：**
<<<代码块_3>>>

---

##### `get_tool_info(tool_name)`
获取有关特定工具的详细信息。

**参数：**
- `tool_name` (str)：工具名称

**返回：** 包含工具元数据、参数和文档的字典

**示例：**
<<<代码块_4>>>

---

## 内置发现工具

这些是帮助寻找生态系统中其他工具的特殊工具。

### 工具_查找器

基于嵌入的语义搜索工具。需要 GPU。

<<<代码块_5>>>

**参数：**
- `description` (str)：所需功能的自然语言描述
- `limit` (int): 返回的最大工具数量

**返回：** 具有相似度分数的相关工具列表

---

### 工具_Finder_LLM

基于 LLM 的语义搜索工具。无需 GPU。

<<<代码块_6>>>

**参数：**
- `description` (str): 自然语言查询
- `limit` (int): 返回的工具的最大数量

**返回：**相关工具列表

---

### 工具_查找器_关键字

通过工具名称和描述进行基于关键字的快速搜索。

```python
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {
        "description": "pathway enrichment",
        "limit": 10
    }
})
```

**参数：**
- `description` (str)：要搜索的关键字
- `limit` (int): 返回的最大工具数量

**返回：** 匹配工具列表

---

## 工具输出挂钩

工具结果的后处理挂钩。

### 总结钩子
```python
result = tu.run({
    "name": "some_tool",
    "arguments": {"param": "value"}
},
hooks={
    "summarize": {
        "format": "brief"  # or "detailed"
    }
})
```

### 文件保存挂钩
```python
result = tu.run({
    "name": "some_tool",
    "arguments": {"param": "value"}
},
hooks={
    "save_to_file": {
        "filename": "output.json",
        "format": "json"  # or "csv", "txt"
    }
})
```

---

## 模型上下文协议（MCP）

### 启动 MCP 服务器

命令行界面：
```bash
tooluniverse-smcp
```

这将启动一个 MCP 服务器，该服务器通过模型上下文协议公开所有 ToolUniverse 工具。

**配置：**
- 默认端口：自动分配
- 协议：MCP标准
- 身份验证：本地使用无需身份验证

---

## 集成模块

### OpenRouter 集成

通过 OpenRouter API 访问 100 多个法学硕士：

```python
from tooluniverse import OpenRouterClient

client = OpenRouterClient(api_key="your_key")
response = client.chat("Analyze this protein sequence", model="anthropic/claude-3-5-sonnet")
```

---

## AI-工具交互协议

ToolUniverse 使用标准化协议进行 LLM 工具通信：

**请求格式：**
```json
{
  "name": "tool_name",
  "arguments": {
    "param1": "value1",
    "param2": "value2"
  }
}
```

**回复格式：**
```json
{
  "status": "success",
  "data": { ... },
  "metadata": {
    "execution_time": 1.23,
    "tool_version": "1.0.0"
  }
}
```

---

## 错误处理

```python
try:
    result = tu.run({
        "name": "some_tool",
        "arguments": {"param": "value"}
    })
except ToolNotFoundError as e:
    print(f"Tool not found: {e}")
except InvalidArgumentError as e:
    print(f"Invalid arguments: {e}")
except ToolExecutionError as e:
    print(f"Execution failed: {e}")
```

---

## 类型提示

```python
from typing import Dict, List, Any, Optional

def run_tool(
    tu: ToolUniverse,
    tool_name: str,
    arguments: Dict[str, Any]
) -> Any:
    """Execute a tool with type-safe arguments."""
    return tu.run({
        "name": tool_name,
        "arguments": arguments
    })
```

---

## 最佳实践

1. **初始化一次**：创建单个 ToolUniverse 实例并重用它
2. **提前加载工具**：启动时调用`load_tools()`一次
3. **Cache Tool Info**：存储常用的工具信息
4. **错误处理**：始终将工具执行包装在 try- except 块中
5. **类型验证**：执行前验证参数类型
6. **资源管理**：考虑远程 API 的速率限制
7. **日志记录**：为生产环境启用日志记录