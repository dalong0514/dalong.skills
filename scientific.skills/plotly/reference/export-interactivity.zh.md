<!-- 此文件由机器翻译自 export-interactivity.md -->

# 导出和交互

## 静态图片导出

### 安装

静态图像导出需要Kaleido：

```bash
uv pip install kaleido
```

Kaleido v1+ 需要您的系统上安装 Chrome/Chromium。

### 支持的格式

- **光栅**：PNG、JPEG、WebP
- **矢量**：SVG、PDF

### 写入文件

<<<代码块_1>>>

### 转换为字节

<<<代码块_2>>>

### 自定义导出

<<<代码块_3>>>

### 设置导出默认值

<<<代码块_4>>>

### 导出多个图形

<<<代码块_5>>>

## 交互式 HTML 导出

### 基本导出

<<<代码块_6>>>

### 文件大小控制

```python
# Full library embedded (~5MB file)
fig.write_html('chart.html', include_plotlyjs=True)

# CDN reference (~2KB file, requires internet)
fig.write_html('chart.html', include_plotlyjs='cdn')

# Local reference (requires plotly.min.js in same directory)
fig.write_html('chart.html', include_plotlyjs='directory')

# No library (for embedding in existing HTML with Plotly.js)
fig.write_html('chart.html', include_plotlyjs=False)
```

### HTML 配置

```python
fig.write_html(
    'chart.html',
    config={
        'displayModeBar': True,
        'displaylogo': False,
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'custom_image',
            'height': 800,
            'width': 1200,
            'scale': 2
        }
    }
)
```

### 嵌入模板

```python
# Get only the div (no full HTML structure)
html_div = fig.to_html(
    full_html=False,
    include_plotlyjs='cdn',
    div_id='my-plot'
)

# Use in Jinja2 template
template = """
<html>
<body>
    <h1>My Dashboard</h1>
    {{ plot_div | safe }}
</body>
</html>
"""
```

## 互动功能

### 内置交互

Plotly 图形自动支持：

- **悬停工具提示** - 悬停时显示数据
- **平移和缩放** - 单击并拖动进行平移，滚动进行缩放
- **框/套索选择** - 选择多个点
- **图例切换** - 单击以隐藏/显示痕迹
- **双击** - 重置轴

### 悬停自定义

```python
# Hover mode
fig.update_layout(
    hovermode='closest'  # 'x', 'y', 'closest', 'x unified', False
)

# Custom hover template
fig.update_traces(
    hovertemplate='<b>%{x}</b><br>' +
                  'Value: %{y:.2f}<br>' +
                  'Extra: %{customdata[0]}<br>' +
                  '<extra></extra>'
)

# Hover data in Plotly Express
fig = px.scatter(
    df, x='x', y='y',
    hover_data={
        'extra_col': True,     # Show column
        'x': ':.2f',           # Format column
        'hidden': False        # Hide column
    },
    hover_name='name_column'   # Bold title
)
```

### 点击事件（Dash/FigureWidget）

对于 Web 应用程序，使用 Dash 或 FigureWidget 进行点击处理：

```python
# With FigureWidget in Jupyter
import plotly.graph_objects as go

fig = go.FigureWidget(data=[go.Scatter(x=[1, 2, 3], y=[4, 5, 6])])

def on_click(trace, points, selector):
    print(f'Clicked on points: {points.point_inds}')

fig.data[0].on_click(on_click)
fig
```

### 缩放和平移

```python
# Disable zoom/pan
fig.update_xaxes(fixedrange=True)
fig.update_yaxes(fixedrange=True)

# Set initial zoom
fig.update_xaxes(range=[0, 10])
fig.update_yaxes(range=[0, 100])

# Constrain zoom
fig.update_xaxes(
    range=[0, 10],
    constrain='domain'
)
```

### Rangeslider（时间序列）

```python
fig = px.line(df, x='date', y='value')

# Add rangeslider
fig.update_xaxes(rangeslider_visible=True)

# Customize rangeslider
fig.update_xaxes(
    rangeslider=dict(
        visible=True,
        thickness=0.05,
        bgcolor='lightgray'
    )
)
```

### 范围选择器按钮

```python
fig.update_xaxes(
    rangeselector=dict(
        buttons=list([
            dict(count=1, label='1m', step='month', stepmode='backward'),
            dict(count=6, label='6m', step='month', stepmode='backward'),
            dict(count=1, label='YTD', step='year', stepmode='todate'),
            dict(count=1, label='1y', step='year', stepmode='backward'),
            dict(step='all', label='All')
        ]),
        x=0.0,
        y=1.0,
        xanchor='left',
        yanchor='top'
    )
)
```

### 按钮和下拉菜单

```python
fig.update_layout(
    updatemenus=[
        dict(
            type='buttons',
            direction='left',
            buttons=list([
                dict(
                    args=[{'type': 'scatter'}],
                    label='Scatter',
                    method='restyle'
                ),
                dict(
                    args=[{'type': 'bar'}],
                    label='Bar',
                    method='restyle'
                )
            ]),
            x=0.1,
            y=1.15
        )
    ]
)
```

### 滑块

```python
fig.update_layout(
    sliders=[
        dict(
            active=0,
            steps=[
                dict(
                    method='update',
                    args=[{'visible': [True, False]},
                          {'title': 'Dataset 1'}],
                    label='Dataset 1'
                ),
                dict(
                    method='update',
                    args=[{'visible': [False, True]},
                          {'title': 'Dataset 2'}],
                    label='Dataset 2'
                )
            ],
            x=0.1,
            y=0,
            len=0.9
        )
    ]
)
```

## 动画

### 使用 Plotly Express

```python
fig = px.scatter(
    df, x='gdp', y='life_exp',
    animation_frame='year',     # Animate over this column
    animation_group='country',  # Group animated elements
    size='population',
    color='continent',
    hover_name='country',
    log_x=True,
    range_x=[100, 100000],
    range_y=[25, 90]
)

# Customize animation speed
fig.layout.updatemenus[0].buttons[0].args[1]['frame']['duration'] = 1000
fig.layout.updatemenus[0].buttons[0].args[1]['transition']['duration'] = 500
```

### 使用图形对象

```python
import plotly.graph_objects as go

fig = go.Figure(
    data=[go.Scatter(x=[1, 2], y=[1, 2])],
    layout=go.Layout(
        updatemenus=[dict(
            type='buttons',
            buttons=[dict(label='Play',
                         method='animate',
                         args=[None])]
        )]
    ),
    frames=[
        go.Frame(data=[go.Scatter(x=[1, 2], y=[1, 2])]),
        go.Frame(data=[go.Scatter(x=[1, 2], y=[2, 3])]),
        go.Frame(data=[go.Scatter(x=[1, 2], y=[3, 4])])
    ]
)
```

## 显示数字

### 在 Jupyter 中

```python
# Default renderer
fig.show()

# Specific renderer
fig.show(renderer='notebook')  # or 'jupyterlab', 'colab', 'kaggle'
```

### 在网络浏览器中

```python
fig.show()  # Opens in default browser
```

### 在达世币应用程序中

```python
import dash
from dash import dcc, html
import plotly.express as px

app = dash.Dash(__name__)

fig = px.scatter(df, x='x', y='y')

app.layout = html.Div([
    dcc.Graph(figure=fig)
])

app.run_server(debug=True)
```

### 保存和加载

```python
# Save as JSON
fig.write_json('figure.json')

# Load from JSON
import plotly.io as pio
fig = pio.read_json('figure.json')

# Save as HTML
fig.write_html('figure.html')
```

## 配置选项

### 显示配置

```python
config = {
    'displayModeBar': True,      # Show toolbar
    'displaylogo': False,        # Hide Plotly logo
    'modeBarButtonsToRemove': ['pan2d', 'lasso2d'],  # Remove buttons
    'toImageButtonOptions': {
        'format': 'png',
        'filename': 'custom_image',
        'height': 500,
        'width': 700,
        'scale': 1
    },
    'scrollZoom': True,          # Enable scroll zoom
    'editable': True,            # Enable editing
    'responsive': True           # Responsive sizing
}

fig.show(config=config)
fig.write_html('chart.html', config=config)
```

### 可用的配置选项

- `displayModeBar`：显示/隐藏工具栏（'悬停'、True、False）
- `displaylogo`：显示 Plotly 徽标
- `modeBarButtonsToRemove`：要隐藏的按钮列表
- `modeBarButtonsToAdd`：自定义按钮
- `scrollZoom`：启用滚动缩放
- `doubleClick`：双击行为（'重置'、'自动调整大小'、'重置+自动调整大小'、False）
- `showAxisDragHandles`：显示轴拖动手柄
- `editable`：允许编辑
- `responsive`：响应式调整大小