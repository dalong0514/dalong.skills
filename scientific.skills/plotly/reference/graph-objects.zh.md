<!-- 此文件由机器翻译自 graph-objects.md -->

# 图形对象 - 低级 API

`plotly.graph_objects` 模块通过表示 Plotly 组件的 Python 类提供对图形构造的细粒度控制。

## 核心课程

- **`go.Figure`** - 主要人物容器
- **`go.FigureWidget`** - Jupyter 兼容的交互式小部件
- **轨迹类型** - 40 多种图表类型（散点图、条形图、热图等）
- **布局组件** - 轴、注释、形状等。

## 主要优势

1. **数据验证** - 针对无效属性的有用错误消息
2. **内置文档** - 可通过文档字符串访问
3. **灵活的语法** - 字典或属性访问
4. **便捷方法** - `.add_trace()`、`.update_layout()`等。
5. **神奇下划线符号** - 紧凑的嵌套属性访问
6. **集成 I/O** - `.show()`、`.write_html()`、`.write_image()`

## 基本图形构建

### 创建空图

```python
import plotly.graph_objects as go

fig = go.Figure()
```

### 添加痕迹

<<<代码块_1>>>

## 常见跟踪类型

### 散布（线条和标记）

<<<代码块_2>>>

### 酒吧

<<<代码块_3>>>

### 热图

<<<代码块_4>>>

### 3D 散点图

<<<代码块_5>>>

## 布局配置

### 更新布局

<<<代码块_6>>>

### 神奇下划线表示法

设置嵌套属性的紧凑方法：

```python
# Instead of:
fig.update_layout(title=dict(text='Title', font=dict(size=20)))

# Use underscores:
fig.update_layout(
    title_text='Title',
    title_font_size=20
)
```

### 轴配置

```python
fig.update_xaxes(
    title='X Axis',
    range=[0, 10],
    showgrid=True,
    gridwidth=1,
    gridcolor='lightgray',
    type='log',  # 'linear', 'log', 'date', 'category'
    tickformat='.2f',
    dtick=1  # Tick spacing
)

fig.update_yaxes(
    title='Y Axis',
    zeroline=True,
    zerolinewidth=2,
    zerolinecolor='black'
)
```

## 更新轨迹

```python
# Update all traces
fig.update_traces(
    marker=dict(size=10, opacity=0.7)
)

# Update specific trace
fig.update_traces(
    marker=dict(color='red'),
    selector=dict(name='Line 1')
)

# Update by position
fig.data[0].marker.size = 15
```

## 添加注释

```python
fig.add_annotation(
    x=2, y=5,
    text='Important Point',
    showarrow=True,
    arrowhead=2,
    arrowsize=1,
    arrowwidth=2,
    arrowcolor='red',
    ax=40,  # Arrow x offset
    ay=-40  # Arrow y offset
)
```

## 添加形状

```python
# Rectangle
fig.add_shape(
    type='rect',
    x0=1, y0=2, x1=3, y1=4,
    line=dict(color='red', width=2),
    fillcolor='lightblue',
    opacity=0.3
)

# Line
fig.add_shape(
    type='line',
    x0=0, y0=0, x1=5, y1=5,
    line=dict(color='green', width=2, dash='dash')
)

# Convenience methods for horizontal/vertical lines
fig.add_hline(y=5, line_dash='dash', line_color='red')
fig.add_vline(x=3, line_dash='dot', line_color='blue')
```

## 图结构

数字遵循树形层次结构：

```python
fig = go.Figure(data=[trace1, trace2], layout=go.Layout(...))

# Access via dictionary syntax
fig['layout']['title'] = 'New Title'
fig['data'][0]['marker']['color'] = 'red'

# Or attribute syntax
fig.layout.title = 'New Title'
fig.data[0].marker.color = 'red'
```

## 复杂图表类型

### 烛台

```python
fig.add_trace(go.Candlestick(
    x=df['date'],
    open=df['open'],
    high=df['high'],
    low=df['low'],
    close=df['close'],
    name='Stock Price'
))
```

### 桑基图

```python
fig = go.Figure(data=[go.Sankey(
    node=dict(
        label=['A', 'B', 'C', 'D'],
        color='blue'
    ),
    link=dict(
        source=[0, 1, 0, 2],
        target=[2, 3, 3, 3],
        value=[8, 4, 2, 8]
    )
)])
```

### 表面（3D）

```python
fig = go.Figure(data=[go.Surface(
    z=z_data,  # 2D array
    x=x_data,
    y=y_data,
    colorscale='Viridis'
)])
```

## 使用 DataFrame

从 pandas DataFrame 构建跟踪：

```python
import pandas as pd

df = pd.DataFrame({
    'x': [1, 2, 3, 4],
    'y': [10, 11, 12, 13]
})

fig = go.Figure()
for group_name, group_df in df.groupby('category'):
    fig.add_trace(go.Scatter(
        x=group_df['x'],
        y=group_df['y'],
        name=group_name,
        mode='lines+markers'
    ))
```

## 何时使用图形对象

在以下情况下使用 graph_objects：
- 创建 Plotly Express 中不可用的图表类型
- 从头开始构建复杂的多轨迹图形
- 需要精确控制各个组件
- 创建专门的可视化（3D 网格、等值面、自定义形状）
- 使用混合图表类型构建子图

在以下情况下使用 Plotly Express：
- 快速创建标准图表
- 使用整洁的 DataFrame 数据
- 想要自动样式和图例