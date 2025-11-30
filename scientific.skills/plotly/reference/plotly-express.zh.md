<!-- 此文件由机器翻译自 plotly-express.md -->

# Plotly Express - 高级 API

Plotly Express (px) 是一个高级界面，用于使用最少的代码（通常为 1-5 行）创建数据可视化。

## 安装

```bash
uv pip install plotly
```

## 主要优势

- 常见图表类型的简洁语法
- 自动颜色编码和图例
- 与 pandas DataFrames 无缝协作
- 布局和样式的智能默认值
- 返回 graph_objects.Figure 以进行进一步定制

## 基本使用模式

<<<代码块_1>>>

## 40+ 图表类型

### 基本图表
- `px.scatter()` - 带有可选趋势线的散点图
- `px.line()` - 时间序列折线图
- `px.bar()` - 条形图（垂直/水平）
- `px.area()` - 面积图
- `px.pie()` - 饼图

### 统计图表
- `px.histogram()` - 带自动分箱的直方图
- `px.box()` - 分布的箱线图
- `px.violin()` - 小提琴图
- `px.strip()` - 带状图
- `px.ecdf()` - 经验累积分布

### 地图
- `px.scatter_geo()` - 地理散点图
- `px.choropleth()` - 等值线地图
- `px.scatter_mapbox()` - Mapbox 散点图
- `px.density_mapbox()` - 地图上的密度热图

### 专业化
- `px.sunburst()` - 层次旭日图
- `px.treemap()` - 树形图可视化
- `px.funnel()` - 漏斗图
- `px.parallel_coordinates()` - 平行坐标
- `px.scatter_matrix()` - 散布矩阵 (SPLOM)
- `px.density_heatmap()` - 2D 密度热图
- `px.density_contour()` - 密度等值线

### 3D 图表
- `px.scatter_3d()` - 3D 散点图
- `px.line_3d()` - 3D 线图

## 常用参数

所有 Plotly Express 函数都支持这些样式参数：

<<<代码块_2>>>

## 数据格式

Plotly Express 适用于：
- **长格式数据**（整洁）：每个观察一行
- **宽格式数据**：多列作为单独的迹线

<<<代码块_3>>>

## 趋势线

添加统计趋势线到散点图：

<<<代码块_4>>>

## 分面（子图）

自动创建多面图：

<<<代码块_5>>>

## 动画

创建动画可视化：

<<<代码块_6>>>

## 悬停数据

自定义悬停工具提示：

```python
fig = px.scatter(
    df, x="x", y="y",
    hover_data={
        "extra_col": True,      # Add column
        "x": ":.2f",            # Format existing
        "hidden_col": False     # Hide column
    },
    hover_name="name_column"    # Bold title in hover
)
```

## 进一步定制

Plotly Express 返回一个可以进一步自定义的 `graph_objects.Figure`：

```python
fig = px.scatter(df, x="x", y="y")

# Use graph_objects methods
fig.update_layout(
    title="Custom Title",
    xaxis_title="X Axis",
    font=dict(size=14)
)

fig.update_traces(
    marker=dict(size=10, opacity=0.7)
)

fig.add_hline(y=0, line_dash="dash")
```

## 何时使用 Plotly Express

在以下情况下使用 Plotly Express：
- 快速创建标准图表类型
- 使用 pandas DataFrames
- 需要自动颜色/尺寸编码
- 想要使用最少的代码进行合理的默认设置

在以下情况下使用 graph_objects：
- 构建不以 px 为单位的自定义图表类型
- 需要对每个元素进行细粒度的控制
- 创建复杂的多轨迹图形
- 构建专业的可视化