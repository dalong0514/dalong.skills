<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：报告实验室
描述：“PDF 生成工具包。创建发票、报告、证书、表格、图表、表格、条形码、QR 码、Canvas/Platypus API，实现专业文档自动化。”
---

# ReportLab PDF 生成

## 概述

ReportLab 是一个功能强大的 Python 库，用于以编程方式生成 PDF。使用表格、图表、图像和交互式表单创建从简单文档到复杂报告的任何内容。

**两种主要方法：**
- **Canvas API**（低级）：使用基于坐标的定位直接绘图 - 用于精确布局
- **Platypus**（高级）：具有自动分页符的流畅文档布局 - 用于多页文档

**核心能力：**
- 具有丰富格式和自定义字体的文本
- 具有复杂样式和单元格跨度的表格
- 图表（条形图、折线图、饼图、面积图、散点图）
- 条形码和 QR 码（Code128、EAN、QR 等）
- 具有透明度的图像
- PDF 功能（链接、书签、表格、加密）

## 选择正确的方法

### 在以下情况下使用 Canvas API：
- 创建标签、名片、证书
- 精确定位至关重要（x、y 坐标）
- 单页文档或简单布局
- 绘制图形、形状和定制设计
- 在特定位置添加条形码或二维码

### 在以下情况下使用 Platypus：
- 创建多页文档（报告、文章、书籍）
- 内容应自动跨页面流动
- 需要在每页上重复的页眉/页脚
- 处理可以跨页拆分的段落
- 使用目录构建复杂的文档

### 在以下情况下同时使用：
- 复杂的报告既需要流畅的内容又需要精确的定位
- 将页眉/页脚添加到 Platypus 文档（使用 Canvas 的 `onPage` 回调）
- 在流动文档（Platypus）中嵌入自定义图形（Canvas）

## 快速入门示例

### 简单画布文档

```python
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

c = canvas.Canvas("output.pdf", pagesize=letter)
width, height = letter

# Draw text
c.setFont("Helvetica-Bold", 24)
c.drawString(inch, height - inch, "Hello ReportLab!")

# Draw a rectangle
c.setFillColorRGB(0.2, 0.4, 0.8)
c.rect(inch, 5*inch, 4*inch, 2*inch, fill=1)

# Save
c.showPage()
c.save()
```

### 简单的鸭嘴兽文档

<<<代码块_1>>>

## 常见任务

### 创建表

表格可与 Canvas（通过 Drawing）和 Platypus（作为 Flowable）配合使用：

<<<代码块_2>>>

**详细表格参考：** 请参阅 `references/tables_reference.md` 了解单元格跨度、边框、对齐方式和高级样式。

### 创建图表

图表使用图形框架，可以添加到 Canvas 和 Platypus 中：

<<<代码块_3>>>

**可用图表类型：** 条形图（垂直/水平）、折线图、饼图、面积图、散点图
**详细图表参考：** 请参阅 `references/charts_reference.md` 了解所有图表类型、样式、图例和自定义。

### 添加条形码和二维码

<<<代码块_4>>>

**支持的格式：** Code128、Code39、EAN-13、EAN-8、UPC-A、ISBN、QR、Data Matrix 等 20 多种
**详细条形码参考：** 请参阅 `references/barcodes_reference.md` 了解所有格式和使用示例。

### 使用文本和字体

<<<代码块_5>>>

**使用自定义字体：**

<<<代码块_6>>>

**详细文本参考：** 有关段落样式、字体系列、亚洲语言、希腊字母和格式，请参阅 `references/text_and_fonts.md`。

### 添加图像

```python
from reportlab.platypus import Image
from reportlab.lib.units import inch

# In Platypus
img = Image('photo.jpg', width=4*inch, height=3*inch)
story.append(img)

# Maintain aspect ratio
img = Image('photo.jpg', width=4*inch, height=3*inch, kind='proportional')

# In Canvas
c.drawImage('photo.jpg', x, y, width=4*inch, height=3*inch)

# With transparency (mask white background)
c.drawImage('logo.png', x, y, mask=[255,255,255,255,255,255])
```

### 创建表格

```python
from reportlab.pdfgen import canvas
from reportlab.lib.colors import black, white, lightgrey

c = canvas.Canvas("form.pdf")

# Text field
c.acroForm.textfield(
    name="name",
    tooltip="Enter your name",
    x=100, y=700,
    width=200, height=20,
    borderColor=black,
    fillColor=lightgrey,
    forceBorder=True
)

# Checkbox
c.acroForm.checkbox(
    name="agree",
    x=100, y=650,
    size=20,
    buttonStyle='check',
    checked=False
)

# Dropdown
c.acroForm.choice(
    name="country",
    x=100, y=600,
    width=150, height=20,
    options=[("United States", "US"), ("Canada", "CA")],
    forceBorder=True
)

c.save()
```

**详细 PDF 功能参考：** 请参阅 `references/pdf_features.md` 了解表单、链接、书签、加密和元数据。

### 页眉和页脚

对于 Platypus 文档，使用页面回调：

```python
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame

def add_header_footer(canvas, doc):
    """Called on each page"""
    canvas.saveState()

    # Header
    canvas.setFont('Helvetica', 9)
    canvas.drawString(inch, height - 0.5*inch, "Document Title")

    # Footer
    canvas.drawRightString(width - inch, 0.5*inch, f"Page {doc.page}")

    canvas.restoreState()

# Set up document
doc = BaseDocTemplate("output.pdf")
frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
template = PageTemplate(id='normal', frames=[frame], onPage=add_header_footer)
doc.addPageTemplates([template])

# Build with story
doc.build(story)
```

## 帮助脚本

此技能包括用于常见任务的帮助程序脚本：

### 快速文档生成器

使用 `scripts/quick_document.py` 快速创建文档：

```python
from scripts.quick_document import create_simple_document, create_styled_table

# Simple document from content blocks
content = [
    {'type': 'heading', 'content': 'Introduction'},
    {'type': 'paragraph', 'content': 'Your text here...'},
    {'type': 'bullet', 'content': 'Bullet point'},
]

create_simple_document("output.pdf", "My Document", content_blocks=content)

# Styled tables with presets
data = [['Header1', 'Header2'], ['Data1', 'Data2']]
table = create_styled_table(data, style_name='striped')  # 'default', 'striped', 'minimal', 'report'
```

## 模板示例

`assets/` 中的完整工作示例：

### 发票模板

`assets/invoice_template.py` - 专业发票包含：
- 公司和客户信息
- 包含计算的详细表格
- 税金和总额
- 条款和注释
- 标志放置

```python
from assets.invoice_template import create_invoice

create_invoice(
    filename="invoice.pdf",
    invoice_number="INV-2024-001",
    invoice_date="January 15, 2024",
    due_date="February 15, 2024",
    company_info={'name': 'Acme Corp', 'address': '...', 'phone': '...', 'email': '...'},
    client_info={'name': 'Client Name', ...},
    items=[
        {'description': 'Service', 'quantity': 1, 'unit_price': 500.00},
        ...
    ],
    tax_rate=0.08,
    notes="Thank you for your business!",
)
```

### 报告模板

`assets/report_template.py` - 多页业务报告，其中包含：
- 封面页
- 目录
- 多个部分与小节
- 图表和表格
- 页眉和页脚

```python
from assets.report_template import create_report

report_data = {
    'title': 'Quarterly Report',
    'subtitle': 'Q4 2023',
    'author': 'Analytics Team',
    'sections': [
        {
            'title': 'Executive Summary',
            'content': 'Report content...',
            'table_data': {...},
            'chart_data': {...}
        },
        ...
    ]
}

create_report("report.pdf", report_data)
```

## 参考文档

按功能组织的综合 API 参考：

- **`references/canvas_api.md`** - 低级画布：绘制基元、坐标、变换、状态管理、图像、路径
- **`references/platypus_guide.md`** - 高级 Platypus：文档模板、框架、流程、页面布局、目录
- **`references/text_and_fonts.md`** - 文本格式：段落样式、内联标记、自定义字体、亚洲语言、项目符号、序列
- **`references/tables_reference.md`** - 表格：创建、样式、单元格跨度、边框、对齐方式、颜色、渐变
- **`references/charts_reference.md`** - 图表：所有图表类型、数据处理、轴、图例、颜色、渲染
- **`references/barcodes_reference.md`** - 条形码：Code128、QR 码、EAN、UPC、邮政编码和 20 多种格式
- **`references/pdf_features.md`** - PDF 功能：链接、书签、表单、加密、元数据、页面转换

## 最佳实践

### 坐标系（画布）
- 原点 (0, 0) 是**左下角**（不是左上角）
- Y轴点**向上**
- 单位为**点**（72 点 = 1 英寸）
- 始终明确指定页面大小

```python
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

width, height = letter
margin = inch

# Top of page
y_top = height - margin

# Bottom of page
y_bottom = margin
```

### 选择页面大小

```python
from reportlab.lib.pagesizes import letter, A4, landscape

# US Letter (8.5" x 11")
pagesize=letter

# ISO A4 (210mm x 297mm)
pagesize=A4

# Landscape
pagesize=landscape(letter)

# Custom
pagesize=(6*inch, 9*inch)
```

### 性能提示

1. **使用 `drawImage()` 而不是 `drawInlineImage()`** - 缓存图像以供重用
2. **启用大文件压缩：** `canvas.Canvas("file.pdf", pageCompression=1)`
3. **重用样式** - 创建一次，在整个文档中使用
4. **使用 Forms/XObjects** 进行重复图形

### 常见模式

**在画布上居中文本：**
```python
text = "Centered Text"
text_width = c.stringWidth(text, "Helvetica", 12)
x = (width - text_width) / 2
c.drawString(x, y, text)

# Or use built-in
c.drawCentredString(width/2, y, text)
```

**Platypus 中的分页符：**
```python
from reportlab.platypus import PageBreak

story.append(PageBreak())
```

**将内容放在一起（不拆分）：**
```python
from reportlab.platypus import KeepTogether

story.append(KeepTogether([
    heading,
    paragraph1,
    paragraph2,
]))
```

**替代行颜色：**
```python
style = TableStyle([
    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightgrey]),
])
```

## 故障排除

**文本重叠或消失：**
- 检查 Y 坐标 - 记住原点是左下角
- 确保文本适合页面边界
- 验证`leading`（行距）大于`fontSize`

**表格不适合页面：**
- 减少列宽
- 减小字体大小
- 使用横向方向
- 使用 `repeatRows` 启用表拆分

**条形码未扫描：**
- 增加`barHeight`（最小尝试0.5英寸）
- 为安静区域设置`quiet=1`
- 测试打印质量（建议 300+ DPI）
- 验证条形码类型的数据格式

**未找到字体：**
- 使用 `pdfmetrics.registerFont()` 注册 TrueType 字体
- 使用与注册完全相同的字体系列名称
- 检查字体文件路径是否正确

**图像有白色背景：**
- 使用`mask`参数使白色透明
- 提供遮罩的 RGB 范围：`mask=[255,255,255,255,255,255]`
- 或者使用带有 alpha 通道的 PNG

## 工作流程示例

### 创建发票

1. 从 `assets/invoice_template.py` 中的发票模板开始
2. 自定义公司信息、Logo路径
3. 添加商品的描述、数量、价格
4. 设置税率（如果适用）
5.添加备注和付款条件
6. 生成PDF

### 创建报告

1. 从 `assets/report_template.py` 的报告模板开始
2. 定义标题和内容的部分
3.使用`create_styled_table()`添加数据表
4.使用图形框架添加图表
5. 使用 `doc.multiBuild(story)` 构建 TOC

### 创建证书

1.使用Canvas API进行精准定位
2.加载自定义字体以实现优雅的排版
3.添加边框图形或图像背景
4. 定位文本元素（姓名、日期、成就）
5. 可选：添加二维码进行验证

### 使用条形码创建标签

1.使用具有自定义页面尺寸（标签尺寸）的Canvas
2.计算每页多个标签的网格位置
3.绘制标签内容（文字、图像）
4. 在特定位置添加条形码
5. 在标签或网格之间使用 `showPage()`

## 安装

```bash
uv pip install reportlab

# For image support
uv pip install pillow

# For charts
uv pip install reportlab[renderPM]

# For barcode support (included in reportlab)
# QR codes require: uv pip install qrcode
```

## 何时使用此技能

该技能应该在以下情况下使用：
- 以编程方式生成 PDF 文档
- 创建发票、收据或账单文件
- 使用表格和图表构建报告
- 生成证书、徽章或凭证
- 创建带有条形码的运输标签或产品标签
- 设计表格或可填写的 PDF
- 生成格式一致的多页文档
- 将数据转换为 PDF 格式以供存档或分发
- 创建需要精确定位的自定义布局

该技能为所有 ReportLab 功能提供全面的指导，从简单的文档到带有图表、表格和交互元素的复杂的多页报告。