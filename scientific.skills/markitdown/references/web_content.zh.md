<!-- 此文件由机器翻译自 web_content.md -->

# 网页内容提取参考

本文档提供有关从 HTML、YouTube、EPUB 和其他基于 Web 的格式提取内容的详细信息。

## HTML 转换

将 HTML 文件和网页转换为干净的 Markdown 格式。

### 基本 HTML 转换

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("webpage.html")
print(result.text_content)
```

### HTML 处理功能

**保留的内容：**
- 标题（`<h1>` → `#`、`<h2>` → `##` 等）
- 段落和文本格式
- 链接 (`<a>` → `[text](url)`)
- 列表（有序和无序）
- 表格 → Markdown 表格
- 代码块和内联代码
- 强调（粗体、斜体）

**删除了什么：**
- 脚本和样式
- 导航元素
- 广告内容
- 样板标记
- HTML 注释

### 来自 URL 的 HTML

直接从 URL 转换网页：

<<<代码块_1>>>

### 干净的网络文章提取

从网络文章中提取主要内容：

<<<代码块_2>>>

### HTML 与图像

包含图像的 HTML 文件可以通过 LLM 描述进行增强：

<<<代码块_3>>>

## YouTube 脚本

从 YouTube 视频中提取视频转录。

### 基本 YouTube 转换

<<<代码块_4>>>

### YouTube 安装

<<<代码块_5>>>

这将安装 `youtube-transcript-api` 依赖项。

### YouTube URL 格式

MarkItDown 支持各种 YouTube URL 格式：
- `https://www.youtube.com/watch?v=VIDEO_ID`
- `https://youtu.be/VIDEO_ID`
- `https://www.youtube.com/embed/VIDEO_ID`
- `https://m.youtube.com/watch?v=VIDEO_ID`

### YouTube 脚本功能

**内含内容：**
- 完整的视频文字记录
- 时间戳（可选，取决于可用性）
- 视频元数据（标题、描述）
- 可用语言的字幕

**转录语言：**
<<<代码块_6>>>

### YouTube 播放列表处理

处理播放列表中的多个视频：

```python
from markitdown import MarkItDown

md = MarkItDown()

video_ids = [
    "VIDEO_ID_1",
    "VIDEO_ID_2",
    "VIDEO_ID_3"
]

transcripts = []
for vid_id in video_ids:
    url = f"https://youtube.com/watch?v={vid_id}"
    result = md.convert(url)
    transcripts.append({
        'video_id': vid_id,
        'transcript': result.text_content
    })
```

### YouTube 用例

**内容分析：**
- 无需观看即可分析视频内容
- 从教程中提取关键信息
- 建立可搜索的转录数据库

**研究：**
- 处理面试笔录
- 提取讲座内容
- 分析演示内容

**辅助功能：**
- 生成视频内容的文本版本
- 创建可搜索的视频档案

### YouTube 限制

- 要求视频有字幕/文字记录
- 自动生成的字幕可能存在转录错误
- 某些视频可能会禁用转录访问
- 速率限制可能适用于批量处理

## EPUB 书籍

将 EPUB 电子书转换为 Markdown 格式。

### 基本 EPUB 转换

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("book.epub")
print(result.text_content)
```

### EPUB 处理功能

**提取的内容：**
- 书籍文字内容
- 章节结构
- 标题和格式
- 目录
- 脚注和参考文献

**保留的内容：**
- 标题层次结构
- 文本强调（粗体、斜体）
- 链接和参考
- 列表和表格

### 带图像的 EPUB

EPUB 文件通常包含图像（封面、图表、插图）：

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
result = md.convert("illustrated_book.epub")
```

### EPUB 用例

**研究：**
- 将教科书转换为可搜索的格式
- 提取内容进行分析
- 建立数字图书馆

**内容处理：**
- 准备LLM培训数据书籍
- 转换为不同的格式
- 创建摘要和摘录

**辅助功能：**
- 转换为更易于访问的格式
- 为屏幕阅读器提取文本
- 文本转语音的过程

## RSS 源

处理 RSS 源以提取文章内容。

### 基本 RSS 处理

```python
from markitdown import MarkItDown
import feedparser

md = MarkItDown()

# Parse RSS feed
feed = feedparser.parse("https://example.com/feed.xml")

# Convert each entry
for entry in feed.entries:
    # Save entry HTML
    with open("temp.html", "w") as f:
        f.write(entry.summary)

    result = md.convert("temp.html")
    print(f"## {entry.title}\n\n{result.text_content}\n\n")
```

## 组合的 Web 内容工作流程

### 网页抓取管道

```python
from markitdown import MarkItDown
import requests
from bs4 import BeautifulSoup

md = MarkItDown()

def scrape_and_convert(url):
    """Scrape webpage and convert to Markdown."""
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')

    # Extract main content
    main_content = soup.find('article') or soup.find('main')

    if main_content:
        # Save HTML
        with open("temp.html", "w") as f:
            f.write(str(main_content))

        # Convert to Markdown
        result = md.convert("temp.html")
        return result.text_content

    return None

# Use it
markdown = scrape_and_convert("https://example.com/article")
print(markdown)
```

### YouTube 学习内容提取

```python
from markitdown import MarkItDown

md = MarkItDown()

# Course videos
course_videos = [
    ("https://youtube.com/watch?v=ID1", "Lesson 1: Introduction"),
    ("https://youtube.com/watch?v=ID2", "Lesson 2: Basics"),
    ("https://youtube.com/watch?v=ID3", "Lesson 3: Advanced")
]

course_content = []
for url, title in course_videos:
    result = md.convert(url)
    course_content.append(f"# {title}\n\n{result.text_content}")

# Combine into course document
full_course = "\n\n---\n\n".join(course_content)
with open("course_transcript.md", "w") as f:
    f.write(full_course)
```

### 文档抓取

```python
from markitdown import MarkItDown
import requests
from urllib.parse import urljoin, urlparse

md = MarkItDown()

def scrape_documentation(base_url, page_urls):
    """Scrape multiple documentation pages."""
    docs = []

    for page_url in page_urls:
        full_url = urljoin(base_url, page_url)

        # Fetch page
        response = requests.get(full_url)
        with open("temp.html", "wb") as f:
            f.write(response.content)

        # Convert
        result = md.convert("temp.html")
        docs.append({
            'url': full_url,
            'content': result.text_content
        })

    return docs

# Example usage
base = "https://docs.example.com/"
pages = ["intro.html", "getting-started.html", "api.html"]
documentation = scrape_documentation(base, pages)
```

### EPUB 库处理

```python
from markitdown import MarkItDown
import os

md = MarkItDown()

def process_epub_library(library_path, output_path):
    """Convert all EPUB books in a directory."""
    for filename in os.listdir(library_path):
        if filename.endswith('.epub'):
            epub_path = os.path.join(library_path, filename)

            try:
                result = md.convert(epub_path)

                # Save markdown
                output_file = filename.replace('.epub', '.md')
                output_full = os.path.join(output_path, output_file)

                with open(output_full, 'w') as f:
                    f.write(result.text_content)

                print(f"Converted: {filename}")
            except Exception as e:
                print(f"Failed to convert {filename}: {e}")

# Process library
process_epub_library("books", "markdown_books")
```

## 错误处理

### HTML 转换错误

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("webpage.html")
    print(result.text_content)
except FileNotFoundError:
    print("HTML file not found")
except Exception as e:
    print(f"Conversion error: {e}")
```

### YouTube 脚本错误

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("https://youtube.com/watch?v=VIDEO_ID")
    print(result.text_content)
except Exception as e:
    print(f"Failed to get transcript: {e}")
    # Common issues: No transcript available, video unavailable, network error
```

### EPUB 转换错误

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("book.epub")
    print(result.text_content)
except Exception as e:
    print(f"EPUB processing error: {e}")
    # Common issues: Corrupted file, unsupported DRM, invalid format
```

## 最佳实践

### HTML 处理
- 转换前清理 HTML 以获得更好的结果
- 使用可读性库提取主要内容
- 适当处理不同的编码
- 删除不必要的标记

### YouTube 处理
- 在批处理之前检查成绩单可用性
- 妥善处理 API 速率限制
- 存储成绩单以避免重新获取
- 尊重 YouTube 的服务条款

### EPUB 处理
- 无法处理受 DRM 保护的 EPUB
- 大型 EPUB 可能需要更多内存
- 某些格式可能无法完美翻译
- 首先使用代表性样品进行测试

### 网页抓取道德
- 尊重robots.txt
- 在请求之间添加延迟
- 在用户代理中识别您的抓取工具
- 缓存结果以最大程度地减少请求
- 遵守网站服务条款