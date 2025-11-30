<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 情节
描述：Python 的交互式科学和统计数据可视化库。在创建图表、绘图或可视化效果时使用，包括散点图、折线图、条形图、热图、3D 图、地理地图、统计分布、财务图表和仪表板。支持快速可视化（Plotly Express）和细粒度定制（图形对象）。输出交互式 HTML 或静态图像（PNG、PDF、SVG）。
---

# 情节

Python 图形库，用于使用 40 多种图表类型创建交互式、出版质量的可视化效果。

## 快速入门

安装情节：
```bash
uv pip install plotly
```

Plotly Express（高级 API）的基本用法：
<<<代码块_1>>>

## 在 API 之间进行选择

### 使用 Plotly Express (px)
对于具有合理默认值的快速、标准可视化：
- 使用 pandas DataFrames
- 创建常见图表类型（散点图、折线图、条形图、直方图等）
- 需要自动颜色编码和图例
- 想要最少的代码（1-5行）

请参阅 [reference/plotly-express.md](reference/plotly-express.md) 以获取完整指南。

### 使用图形对象 (go)
对于细粒度控制和自定义可视化：
- Plotly Express 中没有的图表类型（3D 网格、等值面、复杂的金融图表）
- 从头开始构建复杂的多轨迹图形
- 需要精确控制各个组件
- 使用自定义形状和注释创建专门的可视化效果

请参阅 [reference/graph-objects.md](reference/graph-objects.md) 以获取完整指南。

**注意：** Plotly Express 返回图形对象Figure，因此您可以组合方法：
<<<代码块_2>>>

## 核心能力

### 1. 图表类型

Plotly 支持按类别组织的 40 多种图表类型：

**基本图表：**散点图、折线图、条形图、饼图、面积图、气泡图

**统计图表：**直方图、箱线图、小提琴图、分布图、误差线

**科学图表：** 热图、等高线、三元、图像显示

**财务图表：** 烛台图、OHLC、瀑布图、漏斗图、时间序列

**地图：** 散点图、等值线图、密度图（地理可视化）

**3D 图表：** scatter3d、表面、网格、圆锥体、体积

**专业化：**旭日图、树状图、桑基图、平行坐标、仪表

所有图表类型的详细示例和使用请参见[reference/chart-types.md](reference/chart-types.md)。

### 2. 布局和样式

**子图：** 创建具有共享轴的多图图形：
<<<代码块_3>>>

**模板：** 应用协调的样式：
<<<代码块_4>>>

**定制：** 控制外观的各个方面：
- 颜色（离散序列、连续尺度）
- 字体和文本
- 轴（范围、刻度、网格）
- 传奇
- 边距和尺寸
- 注释和形状

有关完整的布局和样式选项，请参阅 [reference/layouts-styling.md](reference/layouts-styling.md)。

### 3.互动性

内置交互功能：
- 带有可定制数据的悬停工具提示
- 平移和缩放
- 图例切换
- 框/套索选择
- 时间序列的范围滑块
- 按钮和下拉菜单
- 动画

<<<代码块_5>>>

有关完整的交互指南，请参阅 [reference/export-interactivity.md](reference/export-interactivity.md)。

### 4. 导出选项

**交互式 HTML：**
<<<代码块_6>>>

**静态图像（需要kaleido）：**
```bash
uv pip install kaleido
```

```python
fig.write_image('chart.png')   # PNG
fig.write_image('chart.pdf')   # PDF
fig.write_image('chart.svg')   # SVG
```

有关完整的导出选项，请参阅 [reference/export-interactivity.md](reference/export-interactivity.md)。

## 常见工作流程

### 科学数据可视化

```python
import plotly.express as px

# Scatter plot with trendline
fig = px.scatter(df, x='temperature', y='yield', trendline='ols')

# Heatmap from matrix
fig = px.imshow(correlation_matrix, text_auto=True, color_continuous_scale='RdBu')

# 3D surface plot
import plotly.graph_objects as go
fig = go.Figure(data=[go.Surface(z=z_data, x=x_data, y=y_data)])
```

### 统计分析

```python
# Distribution comparison
fig = px.histogram(df, x='values', color='group', marginal='box', nbins=30)

# Box plot with all points
fig = px.box(df, x='category', y='value', points='all')

# Violin plot
fig = px.violin(df, x='group', y='measurement', box=True)
```

### 时间序列和金融

```python
# Time series with rangeslider
fig = px.line(df, x='date', y='price')
fig.update_xaxes(rangeslider_visible=True)

# Candlestick chart
import plotly.graph_objects as go
fig = go.Figure(data=[go.Candlestick(
    x=df['date'],
    open=df['open'],
    high=df['high'],
    low=df['low'],
    close=df['close']
)])
```

### 多图仪表板

```python
from plotly.subplots import make_subplots
import plotly.graph_objects as go

fig = make_subplots(
    rows=2, cols=2,
    subplot_titles=('Scatter', 'Bar', 'Histogram', 'Box'),
    specs=[[{'type': 'scatter'}, {'type': 'bar'}],
           [{'type': 'histogram'}, {'type': 'box'}]]
)

fig.add_trace(go.Scatter(x=[1, 2, 3], y=[4, 5, 6]), row=1, col=1)
fig.add_trace(go.Bar(x=['A', 'B'], y=[1, 2]), row=1, col=2)
fig.add_trace(go.Histogram(x=data), row=2, col=1)
fig.add_trace(go.Box(y=data), row=2, col=2)

fig.update_layout(height=800, showlegend=False)
```

## 与达世币集成

对于交互式 Web 应用程序，请使用 Dash（Plotly 的 Web 应用程序框架）：

```bash
uv pip install dash
```

```python
import dash
from dash import dcc, html
import plotly.express as px

app = dash.Dash(__name__)

fig = px.scatter(df, x='x', y='y')

app.layout = html.Div([
    html.H1('Dashboard'),
    dcc.Graph(figure=fig)
])

app.run_server(debug=True)
```

## 参考文件

- **[plotly-express.md](reference/plotly-express.md)** - 用于快速可视化的高级 API
- **[graph-objects.md](reference/graph-objects.md)** - 用于细粒度控制的低级 API
- **[chart-types.md](reference/chart-types.md)** - 40 多种图表类型的完整目录及示例
- **[layouts-styling.md](reference/layouts-styling.md)** - 子图、模板、颜色、自定义
- **[export-interactivity.md](reference/export-interactivity.md)** - 导出选项和交互功能

## 其他资源

- 官方文档：https://plotly.com/python/
- API 参考：https://plotly.com/python-api-reference/
- 社区论坛：https://community.plotly.com/