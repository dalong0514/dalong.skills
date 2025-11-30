<!-- 此文件由机器翻译自 text_and_fonts.md -->

# 文本和字体参考

ReportLab 中文本格式、段落样式和字体处理的综合指南。

## 文本编码

**重要提示：** 所有文本输入都应为 UTF-8 编码或 Python Unicode 对象（自 ReportLab 2.0 起）。

```python
# Correct - UTF-8 strings
text = "Hello 世界 مرحبا"
para = Paragraph(text, style)

# For legacy data, convert first
import codecs
decoded_text = codecs.decode(legacy_bytes, 'latin-1')
```

## 段落样式

### 创建样式

<<<代码块_1>>>

### 内置样式

<<<代码块_2>>>

## 段落格式

### 基本段落

<<<代码块_3>>>

### 内联格式化标签

<<<代码块_4>>>

### 字体控制

<<<代码块_5>>>

### 上标和下标

<<<代码块_6>>>

### 希腊字母

```python
text = """
<greek>alpha</greek>, <greek>beta</greek>, <greek>gamma</greek>
<greek>epsilon</greek>, <greek>pi</greek>, <greek>omega</greek>
"""

para = Paragraph(text, normal_style)
```

### 链接

```python
# External link
text = '<link href="https://example.com" color="blue">Click here</link>'

# Internal link (to bookmark)
text = '<link href="#section1" color="blue">Go to Section 1</link>'

# Anchor for internal links
text = '<a name="section1"/>Section 1 Heading'

para = Paragraph(text, normal_style)
```

### 内嵌图像

```python
text = """
Here is an inline image: <img src="icon.png" width="12" height="12" valign="middle"/>
"""

para = Paragraph(text, normal_style)
```

### 换行

```python
text = """
First line<br/>
Second line<br/>
Third line
"""

para = Paragraph(text, normal_style)
```

## 字体处理

### 标准字体

ReportLab 包括 14 种标准 PDF 字体（无需嵌入）：

```python
# Helvetica family
'Helvetica'
'Helvetica-Bold'
'Helvetica-Oblique'
'Helvetica-BoldOblique'

# Times family
'Times-Roman'
'Times-Bold'
'Times-Italic'
'Times-BoldItalic'

# Courier family
'Courier'
'Courier-Bold'
'Courier-Oblique'
'Courier-BoldOblique'

# Symbol and Dingbats
'Symbol'
'ZapfDingbats'
```

### TrueType 字体

```python
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# Register single font
pdfmetrics.registerFont(TTFont('CustomFont', 'CustomFont.ttf'))

# Use in Canvas
canvas.setFont('CustomFont', 12)

# Use in Paragraph style
style = ParagraphStyle('Custom', fontName='CustomFont', fontSize=12)
```

### 字体系列

将相关字体注册为一个系列以获得粗体/斜体支持：

```python
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping

# Register fonts
pdfmetrics.registerFont(TTFont('Vera', 'Vera.ttf'))
pdfmetrics.registerFont(TTFont('VeraBd', 'VeraBd.ttf'))
pdfmetrics.registerFont(TTFont('VeraIt', 'VeraIt.ttf'))
pdfmetrics.registerFont(TTFont('VeraBI', 'VeraBI.ttf'))

# Map family (normal, bold, italic, bold-italic)
addMapping('Vera', 0, 0, 'Vera')       # normal
addMapping('Vera', 1, 0, 'VeraBd')     # bold
addMapping('Vera', 0, 1, 'VeraIt')     # italic
addMapping('Vera', 1, 1, 'VeraBI')     # bold-italic

# Now <b> and <i> tags work with this family
style = ParagraphStyle('VeraStyle', fontName='Vera', fontSize=12)
para = Paragraph("Normal <b>Bold</b> <i>Italic</i> <b><i>Both</i></b>", style)
```

### 字体搜索路径

```python
from reportlab.pdfbase.ttfonts import TTFSearchPath

# Add custom font directory
TTFSearchPath.append('/path/to/fonts/')

# Now fonts in this directory can be found by name
pdfmetrics.registerFont(TTFont('MyFont', 'MyFont.ttf'))
```

### 亚洲语言支持

#### 使用 Adobe 语言包（无嵌入）

```python
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont

# Register CID fonts
pdfmetrics.registerFont(UnicodeCIDFont('HeiseiMin-W3'))    # Japanese
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))    # Chinese (Simplified)
pdfmetrics.registerFont(UnicodeCIDFont('MSung-Light'))     # Chinese (Traditional)
pdfmetrics.registerFont(UnicodeCIDFont('HYSMyeongJo-Medium'))  # Korean

# Use in styles
style = ParagraphStyle('Japanese', fontName='HeiseiMin-W3', fontSize=12)
para = Paragraph("日本語テキスト", style)
```

#### 使用带有亚洲字符的 TrueType 字体

```python
# Register TrueType font with full Unicode support
pdfmetrics.registerFont(TTFont('SimSun', 'simsun.ttc'))

style = ParagraphStyle('Chinese', fontName='SimSun', fontSize=12, wordWrap='CJK')
para = Paragraph("中文文本", style)
```

注意：设置 `wordWrap='CJK'` 以在亚洲语言中正确换行。

## 编号和顺序

使用 `<seq>` 标签自动编号：

```python
# Simple numbering
text = "<seq id='chapter'/> Introduction"  # Outputs: 1 Introduction
text = "<seq id='chapter'/> Methods"       # Outputs: 2 Methods

# Reset counter
text = "<seq id='figure' reset='yes'/>"

# Formatting templates
text = "Figure <seq template='%(chapter)s-%(figure+)s' id='figure'/>"
# Outputs: Figure 1-1, Figure 1-2, etc.

# Multi-level numbering
text = "Section <seq template='%(chapter)s.%(section+)s' id='section'/>"
```

## 项目符号和列表

### 使用项目符号样式

```python
bullet_style = ParagraphStyle(
    'Bullet',
    parent=normal_style,
    leftIndent=20,
    bulletIndent=10,
    bulletText='•',          # Unicode bullet
    bulletFontName='Helvetica',
)

story.append(Paragraph("First item", bullet_style))
story.append(Paragraph("Second item", bullet_style))
story.append(Paragraph("Third item", bullet_style))
```

### 自定义项目符号字符

```python
# Different bullet styles
bulletText='•'     # Filled circle
bulletText='◦'     # Open circle
bulletText='▪'     # Square
bulletText='▸'     # Triangle
bulletText='→'     # Arrow
bulletText='1.'    # Numbers
bulletText='a)'    # Letters
```

## 文本测量

```python
from reportlab.pdfbase.pdfmetrics import stringWidth

# Measure string width
width = stringWidth("Hello World", "Helvetica", 12)

# Check if text fits in available width
max_width = 200
if stringWidth(text, font_name, font_size) > max_width:
    # Text is too wide
    pass
```

## 最佳实践

1. **始终使用 UTF-8** 进行文本输入
2. **设置前导 > fontSize** 以提高可读性（通常为 fontSize + 2）
3. **注册字体系列**以获得适当的粗体/斜体支持
4. **如果显示用户内容，则转义 HTML**：使用 `<` 表示 < 并使用 `>` 表示 >
5. **使用 getSampleStyleSheet()** 作为起点，不要从头开始创建所有样式
6. **如果支持多语言内容，请尽早测试亚洲字体**
7. **中文/日文/韩文文本设置 wordWrap='CJK'**
8. **使用 stringWidth()** 在渲染之前检查文本是否适合
9. **在文档开始时定义一次样式**，并在整个过程中重复使用
10. **为对齐文本启用连字符**：`hyphenationLang='en_US'`（需要 pyphen 包）