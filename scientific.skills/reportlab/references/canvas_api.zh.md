<!-- 此文件由机器翻译自 canvas_api.md -->

# Canvas API 参考

Canvas API 使用基于坐标的绘图对 PDF 生成提供低级、精确的控制。

## 坐标系

- 原点 (0, 0) 位于**左下角**（不像网页图形那样位于左上角）
- X轴指向右，Y轴指向上
- 单位为点（72 点 = 1 英寸）
- 默认页面尺寸为A4；明确指定页面大小以保持一致性

## 基本设置

```python
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.units import inch

# Create canvas
c = canvas.Canvas("output.pdf", pagesize=letter)

# Get page dimensions
width, height = letter

# Draw content
c.drawString(100, 100, "Hello World")

# Finish page and save
c.showPage()  # Complete current page
c.save()      # Write PDF to disk
```

## 文字绘制

### 基本字符串方法
<<<代码块_1>>>

### 文本对象（高级）
对于多行且精确控制的复杂文本操作：

<<<代码块_2>>>

## 绘制基元

### 线路
<<<代码块_3>>>

### 形状
<<<代码块_4>>>

### 贝塞尔曲线
<<<代码块_5>>>

## 路径对象

对于复杂的形状，请使用路径对象：

<<<代码块_6>>>

## 颜色

### RGB（屏幕显示）
```python
from reportlab.lib.colors import red, blue, Color

c.setFillColorRGB(r, g, b)      # r, g, b are 0-1
c.setStrokeColorRGB(r, g, b)
c.setFillColor(red)             # Named colors
c.setStrokeColor(blue)

# Custom with alpha transparency
c.setFillColor(Color(0.5, 0, 0, alpha=0.5))
```

### CMYK（专业印刷）
```python
from reportlab.lib.colors import CMYKColor, PCMYKColor

c.setFillColorCMYK(c, m, y, k)  # 0-1 range
c.setStrokeColorCMYK(c, m, y, k)

# Integer percentages (0-100)
c.setFillColor(PCMYKColor(100, 50, 0, 0))
```

## 线条样式

```python
c.setLineWidth(width)           # Thickness in points
c.setLineCap(mode)              # 0=butt, 1=round, 2=square
c.setLineJoin(mode)             # 0=miter, 1=round, 2=bevel
c.setDash(array, phase)         # e.g., [3, 3] for dotted line
```

## 坐标变换

**重要：** 转换是增量和累积的。

```python
# Translation (move origin)
c.translate(dx, dy)

# Rotation (in degrees, counterclockwise)
c.rotate(theta)

# Scaling
c.scale(xscale, yscale)

# Skewing
c.skew(alpha, beta)
```

### 状态管理
```python
# Save current graphics state
c.saveState()

# ... apply transformations and draw ...

# Restore previous state
c.restoreState()
```

**注意：** 状态无法在 `showPage()` 调用之间保留。

## 图片

```python
from reportlab.lib.utils import ImageReader

# Preferred method (with caching)
c.drawImage(image_source, x, y, width=None, height=None,
            mask=None, preserveAspectRatio=False)

# image_source can be:
# - Filename string
# - PIL Image object
# - ImageReader object

# For transparency, specify RGB mask range
c.drawImage("logo.png", 100, 500, mask=[255, 255, 255, 255, 255, 255])

# Inline (inefficient, no caching)
c.drawInlineImage(image_source, x, y, width=None, height=None)
```

## 页面管理

```python
# Complete current page
c.showPage()

# Set page size for next page
c.setPageSize(size)  # e.g., letter, A4

# Page compression (smaller files, slower generation)
c = canvas.Canvas("output.pdf", pageCompression=1)
```

## 常见模式

### 边距和布局
```python
from reportlab.lib.units import inch
from reportlab.lib.pagesizes import letter

width, height = letter
margin = inch

# Draw within margins
content_width = width - 2*margin
content_height = height - 2*margin

# Text at top margin
c.drawString(margin, height - margin, "Header")

# Text at bottom margin
c.drawString(margin, margin, "Footer")
```

### 页眉和页脚
```python
def draw_header_footer(c, width, height):
    c.saveState()
    c.setFont("Helvetica", 9)
    c.drawString(inch, height - 0.5*inch, "Company Name")
    c.drawRightString(width - inch, 0.5*inch, f"Page {c.getPageNumber()}")
    c.restoreState()

# Call on each page
draw_header_footer(c, width, height)
c.showPage()
```

## 最佳实践

1. **始终指定页面大小** - 不同的平台有不同的默认值
2. **使用变量进行测量** - `margin = inch` 而不是硬编码值
3. **匹配 saveState/restoreState** - 始终平衡这些调用
4. **对工程图进行外部应用变换**以防止线宽缩放
5. **使用drawImage而不是drawInlineImage**以获得重复图像的更好性能
6. **自下而上绘制** - 记住Y轴指向上方