<!-- 此文件由机器翻译自 layouts-styling.md -->

# 布局、样式和定制

## 次要情节

### 创建子图

```python
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# Basic grid
fig = make_subplots(rows=2, cols=2)

# Add traces to specific positions
fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6]), row=1, col=1)
fig.add_trace(go.Bar(x=['A', 'B', 'C'], y=[1, 3, 2]), row=1, col=2)
fig.add_trace(go.Scatter(x=[1, 2, 3], y=[2, 3, 4]), row=2, col=1)
```

### 子图选项

<<<代码块_1>>>

### 混合子图类型

<<<代码块_2>>>

### 自定义子图轴

<<<代码块_3>>>

### 共享色阶

<<<代码块_4>>>

## 模板和主题

### 内置模板

<<<代码块_5>>>

### 自定义模板

<<<代码块_6>>>

## 使用 Plotly Express 进行样式设计

### 内置参数

```python
fig = px.scatter(
    df, x='x', y='y',

    # Dimensions
    width=800,
    height=600,

    # Title
    title='Figure Title',

    # Labels
    labels={'x': 'X Axis Label', 'y': 'Y Axis Label'},

    # Colors
    color='category',
    color_discrete_sequence=px.colors.qualitative.Set2,
    color_discrete_map={'A': 'red', 'B': 'blue'},
    color_continuous_scale='Viridis',

    # Ordering
    category_orders={'category': ['A', 'B', 'C']},

    # Template
    template='plotly_white'
)
```

### 设置默认值

```python
import plotly.express as px

# Session-wide defaults
px.defaults.template = 'plotly_white'
px.defaults.width = 800
px.defaults.height = 600
px.defaults.color_continuous_scale = 'Viridis'
```

## 色阶

### 离散颜色

```python
import plotly.express as px

# Named color sequences
color_sequences = [
    px.colors.qualitative.Plotly,
    px.colors.qualitative.D3,
    px.colors.qualitative.G10,
    px.colors.qualitative.Set1,
    px.colors.qualitative.Pastel,
    px.colors.qualitative.Dark2,
]

fig = px.scatter(df, x='x', y='y', color='category',
                color_discrete_sequence=px.colors.qualitative.Set2)
```

### 连续颜色

```python
# Named continuous scales
continuous_scales = [
    'Viridis', 'Plasma', 'Inferno', 'Magma', 'Cividis',  # Perceptually uniform
    'Blues', 'Greens', 'Reds', 'YlOrRd', 'YlGnBu',       # Sequential
    'RdBu', 'RdYlGn', 'Spectral', 'Picnic',              # Diverging
]

fig = px.scatter(df, x='x', y='y', color='value',
                color_continuous_scale='Viridis')

# Reverse scale
fig = px.scatter(df, x='x', y='y', color='value',
                color_continuous_scale='Viridis_r')

# Custom scale
fig = px.scatter(df, x='x', y='y', color='value',
                color_continuous_scale=['blue', 'white', 'red'])
```

### 颜色条定制

```python
fig.update_coloraxes(
    colorbar=dict(
        title='Value',
        tickmode='linear',
        tick0=0,
        dtick=10,
        len=0.7,           # Length relative to plot
        thickness=20,
        x=1.02             # Position
    )
)
```

## 布局定制

### 标题和字体

```python
fig.update_layout(
    title=dict(
        text='Main Title',
        font=dict(size=24, family='Arial', color='darkblue'),
        x=0.5,              # Center title
        xanchor='center'
    ),

    font=dict(
        family='Arial',
        size=14,
        color='black'
    )
)
```

### 边距和尺寸

```python
fig.update_layout(
    width=1000,
    height=600,

    margin=dict(
        l=50,    # left
        r=50,    # right
        t=100,   # top
        b=50,    # bottom
        pad=10   # padding
    ),

    autosize=True  # Auto-resize to container
)
```

### 背景颜色

```python
fig.update_layout(
    plot_bgcolor='#f0f0f0',   # Plot area
    paper_bgcolor='white'      # Figure background
)
```

### 传奇

```python
fig.update_layout(
    showlegend=True,

    legend=dict(
        title='Legend Title',
        orientation='h',           # 'h' or 'v'
        x=0.5,                     # Position
        y=-0.2,
        xanchor='center',
        yanchor='top',
        bgcolor='rgba(255, 255, 255, 0.8)',
        bordercolor='black',
        borderwidth=1,
        font=dict(size=12)
    )
)
```

### 轴

```python
fig.update_xaxes(
    title='X Axis Title',
    title_font=dict(size=16, family='Arial'),

    # Range
    range=[0, 10],
    autorange=True,  # Auto range

    # Grid
    showgrid=True,
    gridwidth=1,
    gridcolor='lightgray',

    # Ticks
    showticklabels=True,
    tickmode='linear',
    tick0=0,
    dtick=1,
    tickformat='.2f',
    tickangle=-45,

    # Zero line
    zeroline=True,
    zerolinewidth=2,
    zerolinecolor='black',

    # Scale
    type='linear',  # 'linear', 'log', 'date', 'category'
)

fig.update_yaxes(
    title='Y Axis Title',
    # ... same options as xaxes
)
```

### 悬停行为

```python
fig.update_layout(
    hovermode='closest',  # 'x', 'y', 'closest', 'x unified', False
)

# Customize hover template
fig.update_traces(
    hovertemplate='<b>%{x}</b><br>Value: %{y:.2f}<extra></extra>'
)
```

### 注释

```python
fig.add_annotation(
    text='Important Note',
    x=2,
    y=5,
    showarrow=True,
    arrowhead=2,
    arrowsize=1,
    arrowwidth=2,
    arrowcolor='red',
    ax=40,  # Arrow x offset
    ay=-40, # Arrow y offset
    font=dict(size=14, color='black'),
    bgcolor='yellow',
    opacity=0.8
)
```

### 形状

```python
# Rectangle
fig.add_shape(
    type='rect',
    x0=1, y0=2, x1=3, y1=4,
    line=dict(color='red', width=2),
    fillcolor='lightblue',
    opacity=0.3
)

# Circle
fig.add_shape(
    type='circle',
    x0=0, y0=0, x1=1, y1=1,
    line_color='purple'
)

# Convenience methods
fig.add_hline(y=5, line_dash='dash', line_color='red',
              annotation_text='Threshold')
fig.add_vline(x=3, line_dash='dot')
fig.add_vrect(x0=1, x1=2, fillcolor='green', opacity=0.2)
fig.add_hrect(y0=4, y1=6, fillcolor='red', opacity=0.2)
```

## 更新方法

### 更新布局

```python
fig.update_layout(
    title='New Title',
    xaxis_title='X',
    yaxis_title='Y'
)
```

### 更新痕迹

```python
# Update all traces
fig.update_traces(marker=dict(size=10, opacity=0.7))

# Update with selector
fig.update_traces(
    marker=dict(color='red'),
    selector=dict(mode='markers', name='Series 1')
)
```

### 更新轴

```python
fig.update_xaxes(showgrid=True, gridcolor='lightgray')
fig.update_yaxes(type='log')
```

## 响应式设计

```python
# Auto-resize to container
fig.update_layout(autosize=True)

# Responsive in HTML
fig.write_html('plot.html', config={'responsive': True})
```