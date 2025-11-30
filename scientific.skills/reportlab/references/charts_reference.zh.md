<!-- 此文件由机器翻译自 charts_reference.md -->

# 图表和图形参考

在 ReportLab 中创建图表和数据可视化的综合指南。

## 图形架构

ReportLab的图形系统提供独立于平台的绘图：

- **图纸** - 形状和图表的容器
- **形状** - 基元（矩形、圆形、直线、多边形、路径）
- **渲染器** - 转换为 PDF、PostScript、SVG 或位图（PNG、GIF、JPG）
- **坐标系** - Y 轴指向上方（与 PDF 类似，与 Web 图形不同）

## 快速入门

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics import renderPDF

# Create drawing (canvas for chart)
drawing = Drawing(400, 200)

# Create chart
chart = VerticalBarChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 125
chart.data = [[100, 150, 130, 180]]
chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4']

# Add chart to drawing
drawing.add(chart)

# Render to PDF
renderPDF.drawToFile(drawing, 'chart.pdf', 'Chart Title')

# Or add as flowable to Platypus document
story.append(drawing)
```

## 可用的图表类型

### 条形图

<<<代码块_1>>>

### 堆积条形图

<<<代码块_2>>>

### 水平条形图

<<<代码块_3>>>

### 折线图

<<<代码块_4>>>

### 线图（X-Y 图）

<<<代码块_5>>>

### 饼图

<<<代码块_6>>>

### 带侧标签的饼图

```python
from reportlab.graphics.charts.piecharts import Pie

chart = Pie()
# ... set position, data, labels ...

# Side label mode (labels in columns beside pie)
chart.sideLabels = 1
chart.sideLabelsOffset = 0.1  # Distance from pie

# Simple labels (not fancy layout)
chart.simpleLabels = 1
```

### 面积图

```python
from reportlab.graphics.charts.areacharts import HorizontalAreaChart

chart = HorizontalAreaChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Areas stack on top of each other
chart.data = [
    [100, 150, 130, 180],  # Bottom area
    [50, 70, 60, 90],      # Top area
]

chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4']

# Area colors
chart.strands[0].fillColor = colors.lightblue
chart.strands[1].fillColor = colors.pink
```

### 散点图

```python
from reportlab.graphics.charts.lineplots import ScatterPlot

chart = ScatterPlot()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Data points
chart.data = [
    [(1, 2), (2, 3), (3, 5), (4, 4), (5, 6)],  # Series 1
    [(1, 1), (2, 2), (3, 3), (4, 3), (5, 4)],  # Series 2
]

# Hide lines, show points only
chart.lines[0].strokeColor = None
chart.lines[1].strokeColor = None

# Marker symbols
from reportlab.graphics.widgets.markers import makeMarker
chart.lines[0].symbol = makeMarker('Circle')
chart.lines[1].symbol = makeMarker('Square')
```

## 轴配置

### 类别轴 (XCategoryAxis)

对于分类数据（标签，而不是数字）：

```python
# Access via chart
axis = chart.categoryAxis

# Labels
axis.categoryNames = ['Jan', 'Feb', 'Mar', 'Apr']

# Label angle (for long labels)
axis.labels.angle = 45
axis.labels.dx = 0
axis.labels.dy = -5

# Label formatting
axis.labels.fontSize = 10
axis.labels.fontName = 'Helvetica'

# Visibility
axis.visible = 1
```

### 值轴 (YValueAxis)

对于数值数据：

```python
# Access via chart
axis = chart.valueAxis

# Range
axis.valueMin = 0
axis.valueMax = 200
axis.valueStep = 50  # Tick interval

# Or auto-configure
axis.valueSteps = [0, 50, 100, 150, 200]  # Explicit steps

# Label formatting
axis.labels.fontSize = 10
axis.labelTextFormat = '%d%%'  # Add percentage sign

# Grid lines
axis.strokeWidth = 1
axis.strokeColor = colors.black
```

## 样式和定制

### 颜色

```python
from reportlab.lib import colors

# Named colors
colors.blue, colors.red, colors.green, colors.yellow

# RGB
colors.Color(0.5, 0.5, 0.5)  # Grey

# With alpha
colors.Color(1, 0, 0, alpha=0.5)  # Semi-transparent red

# Hex colors
colors.HexColor('#FF5733')
```

### 线条样式

```python
# For line charts
chart.lines[0].strokeColor = colors.blue
chart.lines[0].strokeWidth = 2
chart.lines[0].strokeDashArray = [2, 2]  # Dashed line
```

### 条形标签

```python
# Show values on bars
chart.barLabels.nudge = 5  # Offset from bar top
chart.barLabels.fontSize = 8
chart.barLabelFormat = '%d'  # Number format

# For negative values
chart.barLabels.dy = -5  # Position below bar
```

## 传奇

图表可以有关联的图例：

```python
from reportlab.graphics.charts.legends import Legend

# Create legend
legend = Legend()
legend.x = 350
legend.y = 150
legend.columnMaximum = 10

# Link to chart (share colors)
legend.colorNamePairs = [
    (chart.bars[0].fillColor, 'Series 1'),
    (chart.bars[1].fillColor, 'Series 2'),
]

# Add to drawing
drawing.add(legend)
```

## 绘制形状

### 基本形状

```python
from reportlab.graphics.shapes import (
    Drawing, Rect, Circle, Ellipse, Line, Polygon, String
)
from reportlab.lib import colors

drawing = Drawing(400, 200)

# Rectangle
rect = Rect(50, 50, 100, 50)
rect.fillColor = colors.blue
rect.strokeColor = colors.black
rect.strokeWidth = 1
drawing.add(rect)

# Circle
circle = Circle(200, 100, 30)
circle.fillColor = colors.red
drawing.add(circle)

# Line
line = Line(50, 150, 350, 150)
line.strokeColor = colors.black
line.strokeWidth = 2
drawing.add(line)

# Text
text = String(50, 175, "Label Text")
text.fontSize = 12
text.fontName = 'Helvetica'
drawing.add(text)
```

### 路径（复杂形状）

```python
from reportlab.graphics.shapes import Path

path = Path()
path.moveTo(50, 50)
path.lineTo(100, 100)
path.curveTo(120, 120, 140, 100, 150, 50)
path.closePath()

path.fillColor = colors.lightblue
path.strokeColor = colors.blue
path.strokeWidth = 2

drawing.add(path)
```

## 渲染选项

### 渲染为 PDF

```python
from reportlab.graphics import renderPDF

# Direct to file
renderPDF.drawToFile(drawing, 'output.pdf', 'Chart Title')

# As flowable in Platypus
story.append(drawing)
```

### 渲染到图像

```python
from reportlab.graphics import renderPM

# PNG
renderPM.drawToFile(drawing, 'chart.png', fmt='PNG')

# GIF
renderPM.drawToFile(drawing, 'chart.gif', fmt='GIF')

# JPG
renderPM.drawToFile(drawing, 'chart.jpg', fmt='JPG')

# With specific DPI
renderPM.drawToFile(drawing, 'chart.png', fmt='PNG', dpi=150)
```

### 渲染为 SVG

```python
from reportlab.graphics import renderSVG

renderSVG.drawToFile(drawing, 'chart.svg')
```

## 高级定制

### 检查属性

```python
# List all properties
print(chart.getProperties())

# Dump properties (for debugging)
chart.dumpProperties()

# Set multiple properties
chart.setProperties({
    'width': 400,
    'height': 200,
    'data': [[100, 150, 130]],
})
```

### 系列自定义颜色

```python
# Define color scheme
from reportlab.lib.colors import PCMYKColor

colors_list = [
    PCMYKColor(100, 67, 0, 23),   # Blue
    PCMYKColor(0, 100, 100, 0),   # Red
    PCMYKColor(66, 13, 0, 22),    # Green
]

# Apply to chart
for i, color in enumerate(colors_list):
    chart.bars[i].fillColor = color
```

## 完整示例

### 销售报告条形图

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.legends import Legend
from reportlab.lib import colors

drawing = Drawing(400, 250)

# Create chart
chart = VerticalBarChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Data
chart.data = [
    [120, 150, 180, 200],  # 2023
    [100, 130, 160, 190],  # 2022
]
chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4']

# Styling
chart.bars[0].fillColor = colors.HexColor('#3498db')
chart.bars[1].fillColor = colors.HexColor('#e74c3c')
chart.valueAxis.valueMin = 0
chart.valueAxis.valueMax = 250
chart.categoryAxis.labels.fontSize = 10
chart.valueAxis.labels.fontSize = 10

# Add legend
legend = Legend()
legend.x = 325
legend.y = 200
legend.columnMaximum = 2
legend.colorNamePairs = [
    (chart.bars[0].fillColor, '2023'),
    (chart.bars[1].fillColor, '2022'),
]

drawing.add(chart)
drawing.add(legend)

# Add to story or save
story.append(drawing)
```

### 多线趋势图

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.linecharts import HorizontalLineChart
from reportlab.lib import colors

drawing = Drawing(400, 250)

chart = HorizontalLineChart()
chart.x = 50
chart.y = 50
chart.width = 320
chart.height = 170

# Data
chart.data = [
    [10, 15, 12, 18, 20, 25],  # Product A
    [8, 10, 14, 16, 18, 22],   # Product B
    [12, 11, 13, 15, 17, 19],  # Product C
]

chart.categoryAxis.categoryNames = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun']

# Line styling
chart.lines[0].strokeColor = colors.blue
chart.lines[0].strokeWidth = 2
chart.lines[1].strokeColor = colors.red
chart.lines[1].strokeWidth = 2
chart.lines[2].strokeColor = colors.green
chart.lines[2].strokeWidth = 2

# Axes
chart.valueAxis.valueMin = 0
chart.valueAxis.valueMax = 30
chart.categoryAxis.labels.angle = 0
chart.categoryAxis.labels.fontSize = 9
chart.valueAxis.labels.fontSize = 9

drawing.add(chart)
story.append(drawing)
```

## 最佳实践

1. **为绘图设置明确的尺寸**以确保尺寸一致
2. **位置图表**有足够的边距（x，y距边缘至少30-50）
3. **在整个文档中使用一致的配色方案**
4. **明确设置 valueMin 和 valueMax** 以获得一致的比例
5. **使用真实数据进行测试**以确保标签贴合且不重叠
6. **为多系列图表添加图例**
7. **角度类别标签** 如果它们很长（45°效果很好）
8. **保持简单** - 更少的数据系列更容易阅读
9. **使用适当的图表类型** - 条形图用于比较，线条用于趋势，饼图用于比例
10. **考虑色盲友好的调色板** - 避免红/绿组合