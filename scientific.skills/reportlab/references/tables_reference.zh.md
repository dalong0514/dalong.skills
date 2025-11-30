<!-- 此文件由机器翻译自 tables_reference.md -->

# 表格参考

在 ReportLab 中创建表格并设计样式的综合指南。

## 基本表创建

```python
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors

# Simple data (list of lists or tuples)
data = [
    ['Header 1', 'Header 2', 'Header 3'],
    ['Row 1, Col 1', 'Row 1, Col 2', 'Row 1, Col 3'],
    ['Row 2, Col 1', 'Row 2, Col 2', 'Row 2, Col 3'],
]

# Create table
table = Table(data)

# Add to story
story.append(table)
```

## 表构造函数

<<<代码块_1>>>

### 列宽

<<<代码块_2>>>

## 单元格内容类型

### 文本和换行符

<<<代码块_3>>>

### 段落对象

<<<代码块_4>>>

### 图片

<<<代码块_5>>>

### 嵌套表

<<<代码块_6>>>

## 表格样式

使用命令列表应用样式：

```python
from reportlab.platypus import TableStyle
from reportlab.lib import colors

style = TableStyle([
    # Command format: ('COMMAND', (startcol, startrow), (endcol, endrow), *args)
    ('GRID', (0, 0), (-1, -1), 1, colors.black),  # Grid over all cells
    ('BACKGROUND', (0, 0), (-1, 0), colors.grey),  # Header background
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),  # Header text color
])

table = Table(data)
table.setStyle(style)
```

### 单元坐标系

- 列和行的索引为 0：`(col, row)`
- 负索引从末尾开始计数：`-1` 是最后一列/行
- `(0, 0)` 是左上角单元格
- `(-1, -1)` 是右下单元格

```python
# Examples:
(0, 0), (2, 0)      # First three cells of header row
(0, 1), (-1, -1)    # All cells except header
(0, 0), (-1, -1)    # Entire table
```

## 样式命令

### 文本格式

```python
style = TableStyle([
    # Font name
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),

    # Font size
    ('FONTSIZE', (0, 0), (-1, -1), 10),

    # Text color
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
    ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),

    # Combined font command
    ('FONT', (0, 0), (-1, 0), 'Helvetica-Bold', 12),  # name, size
])
```

### 对齐

```python
style = TableStyle([
    # Horizontal alignment: LEFT, CENTER, RIGHT, DECIMAL
    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
    ('ALIGN', (0, 1), (0, -1), 'LEFT'),      # First column left
    ('ALIGN', (1, 1), (-1, -1), 'RIGHT'),    # Other columns right

    # Vertical alignment: TOP, MIDDLE, BOTTOM
    ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
    ('VALIGN', (0, 0), (-1, 0), 'BOTTOM'),   # Header bottom-aligned
])
```

### 单元格填充

```python
style = TableStyle([
    # Individual padding
    ('LEFTPADDING', (0, 0), (-1, -1), 12),
    ('RIGHTPADDING', (0, 0), (-1, -1), 12),
    ('TOPPADDING', (0, 0), (-1, -1), 6),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 6),

    # Or set all at once by setting each
])
```

### 背景颜色

```python
style = TableStyle([
    # Solid background
    ('BACKGROUND', (0, 0), (-1, 0), colors.blue),
    ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),

    # Alternating row colors
    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightblue]),

    # Alternating column colors
    ('COLBACKGROUNDS', (0, 0), (-1, -1), [colors.white, colors.lightgrey]),
])
```

### 渐变背景

```python
from reportlab.lib.colors import Color

style = TableStyle([
    # Vertical gradient (top to bottom)
    ('BACKGROUND', (0, 0), (-1, 0), colors.blue),
    ('VERTICALGRADIENT', (0, 0), (-1, 0),
     [colors.blue, colors.lightblue]),

    # Horizontal gradient (left to right)
    ('HORIZONTALGRADIENT', (0, 1), (-1, 1),
     [colors.red, colors.yellow]),
])
```

### 线条和边框

```python
style = TableStyle([
    # Complete grid
    ('GRID', (0, 0), (-1, -1), 1, colors.black),

    # Box/outline only
    ('BOX', (0, 0), (-1, -1), 2, colors.black),
    ('OUTLINE', (0, 0), (-1, -1), 2, colors.black),  # Same as BOX

    # Inner grid only
    ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.grey),

    # Directional lines
    ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),   # Header border
    ('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),   # Header bottom
    ('LINEBEFORE', (0, 0), (0, -1), 1, colors.black),  # Left border
    ('LINEAFTER', (-1, 0), (-1, -1), 1, colors.black), # Right border

    # Thickness and color
    ('LINEABOVE', (0, 1), (-1, 1), 0.5, colors.grey),  # Thin grey line
])
```

### 单元格跨越

```python
data = [
    ['Spanning Header', '', ''],           # Span will merge these
    ['A', 'B', 'C'],
    ['D', 'E', 'F'],
]

style = TableStyle([
    # Span 3 columns in first row
    ('SPAN', (0, 0), (2, 0)),

    # Center the spanning cell
    ('ALIGN', (0, 0), (2, 0), 'CENTER'),
])

table = Table(data)
table.setStyle(style)
```

**重要提示：** 跨越的单元格必须包含空字符串 `''`。

### 高级跨越示例

```python
# Span multiple rows and columns
data = [
    ['A', 'B', 'B', 'C'],
    ['A', 'D', 'E', 'F'],
    ['A', 'G', 'H', 'I'],
]

style = TableStyle([
    # Span rows in column 0
    ('SPAN', (0, 0), (0, 2)),  # Merge A cells vertically

    # Span columns in row 0
    ('SPAN', (1, 0), (2, 0)),  # Merge B cells horizontally

    ('GRID', (0, 0), (-1, -1), 1, colors.black),
])
```

## 特殊命令

### 圆角

```python
table = Table(data, cornerRadii=[5, 5, 5, 5])  # [TL, TR, BL, BR]

# Or in style
style = TableStyle([
    ('ROUNDEDCORNERS', [10, 10, 0, 0]),  # Rounded top corners only
])
```

### 没有分裂

防止表在特定位置分裂：

```python
style = TableStyle([
    # Don't split between rows 0 and 2
    ('NOSPLIT', (0, 0), (-1, 2)),
])
```

### 拆分特定样式

当表格拆分时，仅将样式应用于第一部分或最后一部分：

```python
style = TableStyle([
    # Style for first part after split
    ('LINEBELOW', (0, 'splitfirst'), (-1, 'splitfirst'), 2, colors.red),

    # Style for last part after split
    ('LINEABOVE', (0, 'splitlast'), (-1, 'splitlast'), 2, colors.blue),
])
```

## 重复标题

```python
# Repeat first row on each page
table = Table(data, repeatRows=1)

# Repeat first 2 rows
table = Table(data, repeatRows=2)
```

## 完整示例

### 样式报告表

```python
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.units import inch

data = [
    ['Product', 'Quantity', 'Unit Price', 'Total'],
    ['Widget A', '10', '$5.00', '$50.00'],
    ['Widget B', '5', '$12.00', '$60.00'],
    ['Widget C', '20', '$3.00', '$60.00'],
    ['', '', 'Subtotal:', '$170.00'],
]

table = Table(data, colWidths=[2.5*inch, 1*inch, 1*inch, 1*inch])

style = TableStyle([
    # Header row
    ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
    ('FONTSIZE', (0, 0), (-1, 0), 12),
    ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),

    # Data rows
    ('BACKGROUND', (0, 1), (-1, -2), colors.beige),
    ('GRID', (0, 0), (-1, -2), 0.5, colors.grey),
    ('ALIGN', (1, 1), (-1, -1), 'RIGHT'),
    ('ALIGN', (0, 1), (0, -1), 'LEFT'),

    # Total row
    ('BACKGROUND', (0, -1), (-1, -1), colors.lightgrey),
    ('LINEABOVE', (0, -1), (-1, -1), 2, colors.black),
    ('FONTNAME', (2, -1), (-1, -1), 'Helvetica-Bold'),
])

table.setStyle(style)
```

### 交替行颜色

```python
data = [
    ['Name', 'Age', 'City'],
    ['Alice', '30', 'New York'],
    ['Bob', '25', 'Boston'],
    ['Charlie', '35', 'Chicago'],
    ['Diana', '28', 'Denver'],
]

table = Table(data, colWidths=[2*inch, 1*inch, 1.5*inch])

style = TableStyle([
    # Header
    ('BACKGROUND', (0, 0), (-1, 0), colors.darkslategray),
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),

    # Alternating rows (zebra striping)
    ('ROWBACKGROUNDS', (0, 1), (-1, -1),
     [colors.white, colors.lightgrey]),

    # Borders
    ('BOX', (0, 0), (-1, -1), 2, colors.black),
    ('LINEBELOW', (0, 0), (-1, 0), 2, colors.black),

    # Padding
    ('LEFTPADDING', (0, 0), (-1, -1), 12),
    ('RIGHTPADDING', (0, 0), (-1, -1), 12),
    ('TOPPADDING', (0, 0), (-1, -1), 6),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
])

table.setStyle(style)
```

## 最佳实践

1. **显式设置 colWidths** 以保持布局一致
2. **对带有标题的多页表使用repeatRows**
3. **应用填充**以获得更好的可读性（特别是LEFTPADDING和RIGHTPADDING）
4. **使用 ROWBACKGROUNDS** 来交替颜色，而不是为每行设置样式
5. **将空字符串**放入将要跨越的单元格中
6. **尽早使用实际数据量测试分页符**
7. **在单元格中使用段落对象**来处理复杂格式的文本
8. **将 VALIGN 设置为 MIDDLE** 以获得不同行高的更好外观
9. **保持表格简单** - 复杂的嵌套表格很难维护
10. **使用一致的样式** - 定义一次，应用于所有表