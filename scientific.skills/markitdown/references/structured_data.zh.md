<!-- 此文件由机器翻译自 structured_data.md -->

# 结构化数据处理参考

本文档提供有关将结构化数据格式（CSV、JSON、XML）转换为 Markdown 的详细信息。

## CSV 文件

将 CSV（逗号分隔值）文件转换为 Markdown 表。

### 基本 CSV 转换

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.csv")
print(result.text_content)
```

### CSV 到 Markdown 表

CSV 文件会自动转换为 Markdown 表格格式：

**输入 CSV (`data.csv`):**
<<<代码块_1>>>

**输出降价：**
<<<代码块_2>>>

### CSV 转换功能

**保留的内容：**
- 所有列标题
- 所有数据行
- 单元格值（文本和数字）
- 柱式结构

**格式：**
- 标题加粗（Markdown 表格格式）
- 列对齐
- 保留空单元格
- 特殊字符被转义

### 大型 CSV 文件

对于大型 CSV 文件：

<<<代码块_3>>>

**性能考虑：**
- 非常大的文件可能需要一些时间来处理
- 考虑预览前几行进行测试
- 内存使用量随文件大小变化
- 非常宽的表格可能无法在所有 Markdown 查看器中正常显示

### 带有特殊字符的 CSV

包含特殊字符的 CSV 文件会自动处理：

<<<代码块_4>>>

### CSV 分隔符

支持标准 CSV 分隔符：
- 逗号 (`,`) - 标准
- 分号 (`;`) - 欧洲格式中常见
- 选项卡 (`\t`) - TSV 文件

### 命令行 CSV 转换

<<<代码块_5>>>

## JSON 文件

将 JSON 数据转换为可读的 Markdown 格式。

### 基本 JSON 转换

<<<代码块_6>>>

### JSON 格式

JSON 转换为可读的、结构化的 Markdown 格式：

**输入 JSON (`config.json`):**
```json
{
  "name": "MyApp",
  "version": "1.0.0",
  "dependencies": {
    "library1": "^2.0.0",
    "library2": "^3.1.0"
  },
  "features": ["auth", "api", "database"]
}
```

**输出降价：**
```markdown
## Configuration

**name:** MyApp
**version:** 1.0.0

### dependencies
- **library1:** ^2.0.0
- **library2:** ^3.1.0

### features
- auth
- api
- database
```

### JSON 数组处理

JSON 数组转换为列表或表格：

**对象数组：**
```json
[
  {"id": 1, "name": "Alice", "active": true},
  {"id": 2, "name": "Bob", "active": false}
]
```

**转换为表格：**
```markdown
| id | name  | active |
|----|-------|--------|
| 1  | Alice | true   |
| 2  | Bob   | false  |
```

### 嵌套 JSON 结构

嵌套 JSON 使用适当的缩进和层次结构进行转换：

```python
from markitdown import MarkItDown

md = MarkItDown()

# Handles deeply nested structures
result = md.convert("complex_config.json")
print(result.text_content)
```

### JSON 行 (JSONL)

对于 JSON Lines 格式（每行一个 JSON 对象）：

```python
from markitdown import MarkItDown
import json

md = MarkItDown()

# Read JSONL file
with open("data.jsonl", "r") as f:
    for line in f:
        obj = json.loads(line)

        # Convert to JSON temporarily
        with open("temp.json", "w") as temp:
            json.dump(obj, temp)

        result = md.convert("temp.json")
        print(result.text_content)
        print("\n---\n")
```

### 大型 JSON 文件

对于大型 JSON 文件：

```python
from markitdown import MarkItDown

md = MarkItDown()

# Convert large JSON
result = md.convert("large_data.json")

# Save to file
with open("output.md", "w") as f:
    f.write(result.text_content)
```

## XML 文件

将 XML 文档转换为结构化 Markdown。

### 基本 XML 转换

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.xml")
print(result.text_content)
```

### XML 结构保留

XML 转换为 Markdown 并保持层次结构：

**输入 XML (`book.xml`):**
```xml
<?xml version="1.0"?>
<book>
  <title>Example Book</title>
  <author>John Doe</author>
  <chapters>
    <chapter id="1">
      <title>Introduction</title>
      <content>Chapter 1 content...</content>
    </chapter>
    <chapter id="2">
      <title>Background</title>
      <content>Chapter 2 content...</content>
    </chapter>
  </chapters>
</book>
```

**输出降价：**
```markdown
# book

## title
Example Book

## author
John Doe

## chapters

### chapter (id: 1)
#### title
Introduction

#### content
Chapter 1 content...

### chapter (id: 2)
#### title
Background

#### content
Chapter 2 content...
```

### XML 属性

XML 属性在转换中保留：

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.xml")
# Attributes shown as (attr: value) in headings
```

### XML 命名空间

XML 命名空间的处理：

```python
from markitdown import MarkItDown

md = MarkItDown()

# Handles xmlns and namespaced elements
result = md.convert("namespaced.xml")
```

### XML 用例

**配置文件：**
- 将 XML 配置转换为可读格式
- 记录系统配置
- 比较配置文件

**数据交换：**
- 转换 XML API 响应
- 处理 XML 数据源
- 格式之间的转换

**文件处理：**
- 将 DocBook 转换为 Markdown
- 处理SVG描述
- 提取结构化数据

## 结构化数据工作流程

### CSV 数据分析管道

```python
from markitdown import MarkItDown
import pandas as pd

md = MarkItDown()

# Read CSV for analysis
df = pd.read_csv("data.csv")

# Do analysis
summary = df.describe()

# Convert both to Markdown
original = md.convert("data.csv")

# Save summary as CSV then convert
summary.to_csv("summary.csv")
summary_md = md.convert("summary.csv")

print("## Original Data\n")
print(original.text_content)
print("\n## Statistical Summary\n")
print(summary_md.text_content)
```

### JSON API 文档

```python
from markitdown import MarkItDown
import requests
import json

md = MarkItDown()

# Fetch JSON from API
response = requests.get("https://api.example.com/data")
data = response.json()

# Save as JSON
with open("api_response.json", "w") as f:
    json.dump(data, f, indent=2)

# Convert to Markdown
result = md.convert("api_response.json")

# Create documentation
doc = f"""# API Response Documentation

## Endpoint
GET https://api.example.com/data

## Response
{result.text_content}
"""

with open("api_docs.md", "w") as f:
    f.write(doc)
```

### XML 到 Markdown 文档

```python
from markitdown import MarkItDown

md = MarkItDown()

# Convert XML documentation
xml_files = ["config.xml", "schema.xml", "data.xml"]

for xml_file in xml_files:
    result = md.convert(xml_file)

    output_name = xml_file.replace('.xml', '.md')
    with open(f"docs/{output_name}", "w") as f:
        f.write(result.text_content)
```

### 多格式数据处理

```python
from markitdown import MarkItDown
import os

md = MarkItDown()

def convert_structured_data(directory):
    """Convert all structured data files in directory."""
    extensions = {'.csv', '.json', '.xml'}

    for filename in os.listdir(directory):
        ext = os.path.splitext(filename)[1]

        if ext in extensions:
            input_path = os.path.join(directory, filename)
            result = md.convert(input_path)

            # Save Markdown
            output_name = filename.replace(ext, '.md')
            output_path = os.path.join("markdown", output_name)

            with open(output_path, 'w') as f:
                f.write(result.text_content)

            print(f"Converted: {filename} → {output_name}")

# Process all structured data
convert_structured_data("data")
```

### CSV 到 JSON 到 Markdown

```python
import pandas as pd
from markitdown import MarkItDown
import json

md = MarkItDown()

# Read CSV
df = pd.read_csv("data.csv")

# Convert to JSON
json_data = df.to_dict(orient='records')
with open("temp.json", "w") as f:
    json.dump(json_data, f, indent=2)

# Convert JSON to Markdown
result = md.convert("temp.json")
print(result.text_content)
```

### 数据库导出到 Markdown

```python
from markitdown import MarkItDown
import sqlite3
import csv

md = MarkItDown()

# Export database query to CSV
conn = sqlite3.connect("database.db")
cursor = conn.execute("SELECT * FROM users")

with open("users.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow([description[0] for description in cursor.description])
    writer.writerows(cursor.fetchall())

# Convert to Markdown
result = md.convert("users.csv")
print(result.text_content)
```

## 错误处理

### CSV 错误

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("data.csv")
    print(result.text_content)
except FileNotFoundError:
    print("CSV file not found")
except Exception as e:
    print(f"CSV conversion error: {e}")
    # Common issues: encoding problems, malformed CSV, delimiter issues
```

### JSON 错误

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("data.json")
    print(result.text_content)
except Exception as e:
    print(f"JSON conversion error: {e}")
    # Common issues: invalid JSON syntax, encoding issues
```

### XML 错误

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("data.xml")
    print(result.text_content)
except Exception as e:
    print(f"XML conversion error: {e}")
    # Common issues: malformed XML, encoding problems, namespace issues
```

## 最佳实践

### CSV 处理
- 转换前检查分隔符
- 验证编码（推荐UTF-8）
- 如果需要，可以通过流处理处理大文件
- 预览非常宽的表格的输出

### JSON 处理
- 转换前验证 JSON
- 考虑漂亮地打印复杂的结构
- 适当处理循环引用
- 注意大型阵列性能

### XML 处理
- 首先验证 XML 结构
- 一致地处理命名空间
- 考虑使用 XPath 进行选择性提取
- 注意非常深的嵌套

### 数据质量
- 尽可能在转换前清理数据
- 适当处理缺失值
- 验证特殊字符处理
- 使用代表性样品进行测试

### 性能
- 批量处理大文件
- 对非常大的数据集使用流式传输
- 监控内存使用情况
- 适当时缓存转换结果