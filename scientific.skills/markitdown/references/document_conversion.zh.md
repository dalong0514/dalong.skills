<!-- 此文件由机器翻译自 document_conversion.md -->

# 文档转换参考

本文档提供有关使用 MarkItDown 将 Office 文档和 PDF 转换为 Markdown 的详细信息。

## PDF 文件

PDF 转换从 PDF 文档中提取文本、表格和结构。

### 基本 PDF 转换

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("document.pdf")
print(result.text_content)
```

### PDF 与 Azure 文档智能

对于包含表格、表单和复杂布局的复杂 PDF，请使用 Azure 文档智能来增强提取：

<<<代码块_1>>>

**Azure 文档智能的优势：**
- 卓越的表格提取和重建
- 更好地处理多列布局
- 表单字段识别
- 改进了复杂文档中的文本排序

### PDF 处理注意事项

- 扫描的 PDF 需要 OCR（如果安装了 tesseract，则会自动处理）
- 不支持受密码保护的 PDF
- 大型 PDF 可能需要更长的时间来处理
- 尽可能提取矢量图形和嵌入图像

## Word 文档 (DOCX)

Word 文档转换保留标题、段落、列表、表格和超链接。

### 基本 DOCX 转换

<<<代码块_2>>>

### DOCX 结构保留

MarkItDown 保留：
- **标题** → Markdown 标题（`#`、`##` 等）
- **粗体/斜体** → Markdown 强调 (`**bold**`, `*italic*`)
- **列表** → Markdown 列表（有序和无序）
- **表格** → Markdown 表格
- **超链接** → Markdown 链接 `[text](url)`
- **图像** → 参考描述（可以使用 LLM 进行描述）

### 命令行用法

<<<代码块_3>>>

### 带图像的 DOCX

要生成 Word 文档中图像的描述，请使用 LLM 集成：

<<<代码块_4>>>

## PowerPoint 演示文稿 (PPTX)

PowerPoint 转换从幻灯片中提取文本，同时保留结构。

### 基本 PPTX 转换

<<<代码块_5>>>

### PPTX 结构

MarkItDown 将演示文稿处理为：
- 每张幻灯片都成为一个主要部分
- 幻灯片标题变为标题
- 保留要点
- 表格转换为 Markdown 表格
- 附注（如果有）

### 带有图像描述的 PPTX

演示文稿通常包含重要的视觉信息。使用LLM积分来描述图像：

<<<代码块_6>>>

**自定义演示提示：**
- “描述图表和图表及其关键数据点”
- “解释图表及其关系”
- “总结视觉内容以方便访问”

## Excel 电子表格（XLSX、XLS）

Excel 转换将电子表格数据格式化为 Markdown 表格。

### 基本 XLSX 转换

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.xlsx")
print(result.text_content)
```

### 多页工作簿

对于具有多个工作表的工作簿：
- 每张纸都成为一个单独的部分
- 工作表名称用作标题
- 跳过空表
- 计算公式（显示值，而不是公式）

### XLSX 转换详细信息

**保留的内容：**
- 单元格值（文本、数字、日期）
- 表结构（行和列）
- 工作表名称
- 单元格格式（粗体标题）

**未保留的内容：**
- 公式（仅计算值）
- 图表和图形（使用 LLM 集成进行描述）
- 单元格颜色和条件格式
- 评论和注释

### 大型电子表格

对于大型电子表格，请考虑：
- 对于具有许多行/列的文件，处理可能会更慢
- 非常宽的表格可能无法在 Markdown 中很好地格式化
- 如果可能的话考虑过滤或预处理数据

### XLS（旧版 Excel）文件

支持旧版 `.xls` 文件，但需要额外的依赖项：

```bash
pip install 'markitdown[xls]'
```

然后正常使用：
```python
md = MarkItDown()
result = md.convert("legacy_data.xls")
```

## 常见文档转换模式

### 批量文档处理

```python
from markitdown import MarkItDown
import os

md = MarkItDown()

# Process all documents in a directory
for filename in os.listdir("documents"):
    if filename.endswith(('.pdf', '.docx', '.pptx', '.xlsx')):
        result = md.convert(f"documents/{filename}")

        # Save to output directory
        output_name = os.path.splitext(filename)[0] + ".md"
        with open(f"markdown/{output_name}", "w") as f:
            f.write(result.text_content)
```

### 具有混合内容的文档

对于包含多种类型内容（文本、表格、图像）的文档：

```python
from markitdown import MarkItDown
from openai import OpenAI

# Use LLM for image descriptions + Azure for complex tables
client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    docintel_endpoint="YOUR-ENDPOINT",
    docintel_key="YOUR-KEY"
)

result = md.convert("complex_report.pdf")
```

### 错误处理

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("document.pdf")
    print(result.text_content)
except Exception as e:
    print(f"Conversion failed: {e}")
    # Handle specific errors (file not found, unsupported format, etc.)
```

## 输出质量提示

**为了获得最佳结果：**
1.使用 Azure 文档智能处理具有复杂表格的 PDF
2.为具有重要视觉内容的文档启用LLM描述
3. 确保源文档结构良好（正确的标题等）
4. 对于扫描文档，确保良好的扫描质量以保证 OCR 准确性
5. 使用样本文档进行测试以验证输出质量

## 性能考虑因素

**转换速度取决于：**
- 文档大小和复杂性
- 图片数量（尤其是法学硕士描述）
- Azure文档智能的使用
- 可用的系统资源

**优化技巧：**
- 如果不需要图像描述，则禁用 LLM 集成
- 对简单文档使用标准提取（不是 Azure）
- 尽可能并行处理大批量
- 考虑对非常大的文档进行流式传输