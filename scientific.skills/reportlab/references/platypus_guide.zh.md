<!-- 此文件由机器翻译自 platypus_guide.md -->

# Platypus 指南 - 高级页面布局

Platypus（“使用脚本的页面布局和版式”）使用最少的代码为复杂、流畅的文档提供高级文档布局。

## 架构概述

Platypus 采用分层设计：

1. **DocTemplates** - 具有页面格式规则的文档容器
2. **PageTemplates** - 不同页面布局的规范
3. **框架** - 内容流动的区域
4. **Flowables** - 内容元素（段落、表格、图像、间隔符）
5. **Canvas** - 底层渲染引擎（通常是隐藏的）

## 快速入门

```python
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch

# Create document
doc = SimpleDocTemplate("output.pdf", pagesize=letter,
                       rightMargin=72, leftMargin=72,
                       topMargin=72, bottomMargin=18)

# Create story (list of flowables)
story = []
styles = getSampleStyleSheet()

# Add content
story.append(Paragraph("Title", styles['Title']))
story.append(Spacer(1, 0.2*inch))
story.append(Paragraph("Body text here", styles['BodyText']))
story.append(PageBreak())

# Build PDF
doc.build(story)
```

## 核心组件

### 文档模板

#### 简单文档模板
标准文档最常见的模板：

<<<代码块_1>>>

#### BaseDocTemplate（高级）
对于具有多个页面布局的复杂文档：

<<<代码块_2>>>

### 框架

框架定义内容流动的区域：

<<<代码块_3>>>

### 页面模板

使用框架和可选功能定义页面布局：

<<<代码块_4>>>

## 流动性

Flowables 是流过框架的内容元素。

### 常见的 Flowable

<<<代码块_5>>>

### 段落可流动
有关段落的详细用法，请参阅 `text_and_fonts.md`。

<<<代码块_6>>>

### 图像可流动

```python
from reportlab.platypus import Image

# Auto-size to fit
img = Image('photo.jpg')

# Fixed size
img = Image('photo.jpg', width=2*inch, height=2*inch)

# Maintain aspect ratio with max width
img = Image('photo.jpg', width=4*inch, height=3*inch,
           kind='proportional')

story.append(img)
```

### 表可流动
有关表的详细用法，请参阅`tables_reference.md`。

```python
from reportlab.platypus import Table

data = [['Header1', 'Header2'],
        ['Row1Col1', 'Row1Col2'],
        ['Row2Col1', 'Row2Col2']]

table = Table(data, colWidths=[2*inch, 2*inch])
story.append(table)
```

## 页面布局

### 单列文档

```python
doc = SimpleDocTemplate("output.pdf", pagesize=letter)
story = []
# Add flowables...
doc.build(story)
```

### 两列布局

```python
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame

doc = BaseDocTemplate("output.pdf", pagesize=letter)
width, height = letter
margin = inch

# Two side-by-side frames
frame1 = Frame(margin, margin, width/2 - 1.5*margin, height - 2*margin, id='col1')
frame2 = Frame(width/2 + 0.5*margin, margin, width/2 - 1.5*margin, height - 2*margin, id='col2')

template = PageTemplate(id='TwoCol', frames=[frame1, frame2])
doc.addPageTemplates([template])

story = []
# Content flows left column first, then right column
# Add flowables...
doc.build(story)
```

### 多页面模板

```python
from reportlab.platypus import NextPageTemplate

# Define templates
cover_template = PageTemplate(id='Cover', frames=[cover_frame])
body_template = PageTemplate(id='Body', frames=[body_frame])

doc.addPageTemplates([cover_template, body_template])

story = []
# Cover page content
story.append(Paragraph("Cover", title_style))
story.append(NextPageTemplate('Body'))  # Switch to body template
story.append(PageBreak())

# Body content
story.append(Paragraph("Chapter 1", heading_style))
# ...

doc.build(story)
```

## 页眉和页脚

页眉和页脚通过 `onPage` 回调函数添加：

```python
def header_footer(canvas, doc):
    """Draw header and footer on each page"""
    canvas.saveState()

    # Header
    canvas.setFont('Helvetica-Bold', 12)
    canvas.drawCentredString(letter[0]/2, letter[1] - 0.5*inch,
                            "Document Title")

    # Footer
    canvas.setFont('Helvetica', 9)
    canvas.drawString(inch, 0.75*inch, "Left Footer")
    canvas.drawRightString(letter[0] - inch, 0.75*inch,
                          f"Page {doc.page}")

    canvas.restoreState()

# Apply to template
template = PageTemplate(id='Normal', frames=[frame], onPage=header_footer)
```

## 目录

```python
from reportlab.platypus import TableOfContents
from reportlab.lib.styles import ParagraphStyle

# Create TOC
toc = TableOfContents()
toc.levelStyles = [
    ParagraphStyle(name='TOC1', fontSize=14, leftIndent=0),
    ParagraphStyle(name='TOC2', fontSize=12, leftIndent=20),
]

story = []
story.append(toc)
story.append(PageBreak())

# Add entries
story.append(Paragraph("Chapter 1<a name='ch1'/>", heading_style))
toc.addEntry(0, "Chapter 1", doc.page, 'ch1')

# Must call build twice for TOC to populate
doc.build(story)
```

## 文档属性

```python
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.units import inch, cm, mm

# Page sizes
letter  # US Letter (8.5" x 11")
A4      # ISO A4 (210mm x 297mm)
landscape(letter)  # Rotate to landscape

# Units
inch    # 72 points
cm      # 28.35 points
mm      # 2.835 points

# Custom page size
custom_size = (6*inch, 9*inch)
```

## 最佳实践

1. **对大多数文档使用 SimpleDocTemplate** - 它处理常见布局
2. **在调用`doc.build(story)`之前完全构建故事列表**
3. **使用 Spacer** 进行垂直间距而不是空段落
4. **使用 KeepTogether 对相关内容进行分组**，以防止尴尬的分页
5. **尽早使用实际内容量测试分页符**
6. **一致地使用样式** - 创建一次样式，在整个文档中重复使用
7. **在开发期间在框架上设置 showBoundary=1** 以可视化布局
8. **页眉/页脚进入 onPage** 回调，而不是故事中
9. **对于长文档**，使用带有多个页面模板的 BaseDocTemplate
10. **构建目录文档两次**以正确填充目录