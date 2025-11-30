<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：降价
描述：将各种文件格式（PDF、Office 文档、图像、音频、网页内容、结构化数据）转换为针对 LLM 处理优化的 Markdown。将文档转换为 Markdown、从 PDF/Office 文件中提取文本、转录音频、对图像执行 OCR、提取 YouTube 转录本或处理批量文件时使用。支持 20 多种格式，包括 DOCX、XLSX、PPTX、PDF、HTML、EPUB、CSV、JSON、带 OCR 的图像以及带转录的音频。
---

# 标记下来

## 概述

MarkItDown 是一个 Python 实用程序，可将各种文件格式转换为 Markdown 格式，并针对大型语言模型和文本分析管道的使用进行了优化。它保留文档结构（标题、列表、表格、超链接），同时生成干净、高效的 Markdown 输出。

## 何时使用此技能

当用户请求时使用此技能：
- 将文档转换为 Markdown 格式
- 从 PDF、Word、PowerPoint 或 Excel 文件中提取文本
- 对图像执行 OCR 以提取文本
- 将音频文件转录为文本
- 提取 YouTube 视频文字记录
- 将 HTML、EPUB 或网页内容处理为 Markdown
- 将结构化数据（CSV、JSON、XML）转换为可读的 Markdown
- 批量转换多个文件或 ZIP 档案
- 为LLM分析或RAG系统准备文件

## 核心能力

### 1. 文档转换

将 Office 文档和 PDF 转换为 Markdown，同时保留结构。

**支持的格式：**
- PDF 文件（具有可选的 Azure 文档智能集成）
- Word文档（DOCX）
- PowerPoint 演示文稿 (PPTX)
- Excel 电子表格（XLSX、XLS）

**基本用法：**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("document.pdf")
print(result.text_content)
```

**命令行：**
<<<代码块_1>>>

有关文档特定功能的详细文档，请参阅 `references/document_conversion.md`。

### 2.媒体处理

使用 OCR 从图像中提取文本并将音频文件转录为文本。

**支持的格式：**
- 具有 EXIF 元数据提取的图像（JPEG、PNG、GIF 等）
- 带有语音转录的音频文件（需要语音识别）

**带有 OCR 的图像：**
<<<代码块_2>>>

**音频转录：**
<<<代码块_3>>>

有关高级媒体处理选项，请参阅`references/media_processing.md`。

### 3.网页内容提取

将基于 Web 的内容和电子书转换为 Markdown。

**支持的格式：**
- HTML 文件和网页
- YouTube 视频文字记录（通过 URL）
- EPUB 书籍
- RSS 源

**YouTube 文字记录：**
<<<代码块_4>>>

有关网页提取的详细信息，请参阅 `references/web_content.md`。

### 4. 结构化数据处理

将结构化数据格式转换为可读的 Markdown 表。

**支持的格式：**
- CSV 文件
- JSON 文件
- XML 文件

**CSV 到 Markdown 表：**
<<<代码块_5>>>

请参阅 `references/structured_data.md` 以了解特定于格式的选项。

### 5. 高级集成

利用人工智能支持的功能提高转换质量。

**Azure 文档智能：**
为了通过更好的表格提取和布局分析来增强 PDF 处理：
<<<代码块_6>>>

**LLM 支持的图像描述：**
使用 GPT-4o 生成详细的图像描述：
```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
result = md.convert("presentation.pptx")  # Images described with LLM
```

有关集成详细信息，请参阅 `references/advanced_integrations.md`。

### 6. 批处理

一次处理多个文件或整个 ZIP 存档。

**ZIP 文件处理：**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("archive.zip")
print(result.text_content)  # All files converted and concatenated
```

**批处理脚本：**
使用提供的批处理脚本进行目录转换：
```bash
python scripts/batch_convert.py /path/to/documents /path/to/output
```

有关实现详细信息，请参阅`scripts/batch_convert.py`。

## 安装

**完整安装（所有功能）：**
```bash
uv pip install 'markitdown[all]'
```

**模块化安装（特定功能）：**
```bash
uv pip install 'markitdown[pdf]'           # PDF support
uv pip install 'markitdown[docx]'          # Word support
uv pip install 'markitdown[pptx]'          # PowerPoint support
uv pip install 'markitdown[xlsx]'          # Excel support
uv pip install 'markitdown[audio]'         # Audio transcription
uv pip install 'markitdown[youtube]'       # YouTube transcripts
```

**要求：**
- Python 3.10 或更高版本

## 输出格式

MarkItDown 生成针对 LLM 消耗进行优化的干净、代币高效的 Markdown：
- 保留标题、列表和表格
- 维护超链接和格式
- 包括相关元数据（EXIF、文档属性）
- 没有创建临时文件（流方法）

## 常见工作流程

**为 RAG 准备文件：**
```python
from markitdown import MarkItDown

md = MarkItDown()

# Convert knowledge base documents
docs = ["manual.pdf", "guide.docx", "faq.html"]
markdown_content = []

for doc in docs:
    result = md.convert(doc)
    markdown_content.append(result.text_content)

# Now ready for embedding and indexing
```

**文档分析管道：**
```bash
# Convert all PDFs in directory
for file in documents/*.pdf; do
    markitdown "$file" -o "markdown/$(basename "$file" .pdf).md"
done
```

## 插件系统

MarkItDown 支持用于自定义转换逻辑的可扩展插件。为了安全起见，插件默认被禁用：

```python
from markitdown import MarkItDown

# Enable plugins if needed
md = MarkItDown(enable_plugins=True)
```

## 资源

该技能包括每种功能的综合参考文档：

- **references/document_conversion.md** - 详细的 PDF、DOCX、PPTX、XLSX 转换选项
- **references/media_processing.md** - 图像 OCR 和音频转录详细信息
- **references/web_content.md** - HTML、YouTube 和 EPUB 提取
- **references/structured_data.md** - CSV、JSON、XML 转换格式
- **references/advanced_integrations.md** - Azure 文档智能和 LLM 集成
- **scripts/batch_convert.py** - 目录的批处理实用程序