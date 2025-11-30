<!-- 此文件由机器翻译自 chart-types.md -->

# Plotly 图表类型

按类别组织的图表类型综合指南。

## 基本图表

### 散点图

```python
import plotly.express as px
fig = px.scatter(df, x='x', y='y', color='category', size='size')

# With trendlines
fig = px.scatter(df, x='x', y='y', trendline='ols')
```

### 折线图

<<<代码块_1>>>

### 条形图

<<<代码块_2>>>

### 饼图

<<<代码块_3>>>

### 面积图

<<<代码块_4>>>

## 统计图表

### 直方图

<<<代码块_5>>>

### 箱线图

<<<代码块_6>>>

### 小提琴情节

```python
fig = px.violin(df, x='category', y='value', color='group', box=True, points='all')
```

### 带/点图

```python
fig = px.strip(df, x='category', y='value', color='group')
```

### 分布图

```python
# Empirical cumulative distribution
fig = px.ecdf(df, x='value', color='group')

# Marginal distribution
fig = px.scatter(df, x='x', y='y', marginal_x='histogram', marginal_y='box')
```

### 误差线

```python
fig = px.scatter(df, x='x', y='y', error_y='error', error_x='x_error')

# Using graph_objects for custom error bars
import plotly.graph_objects as go
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=[1, 2, 3],
    y=[5, 10, 15],
    error_y=dict(
        type='data',
        array=[1, 2, 3],
        visible=True
    )
))
```

## 科学图表

### 热图

```python
# From matrix data
fig = px.imshow(z_matrix, color_continuous_scale='Viridis')

# With graph_objects
fig = go.Figure(data=go.Heatmap(
    z=z_matrix,
    x=x_labels,
    y=y_labels,
    colorscale='RdBu'
))
```

### 等高线图

```python
# 2D contour
fig = px.density_contour(df, x='x', y='y')

# Filled contour
fig = go.Figure(data=go.Contour(
    z=z_matrix,
    contours=dict(
        coloring='heatmap',
        showlabels=True
    )
))
```

### 三元图

```python
fig = px.scatter_ternary(df, a='component_a', b='component_b', c='component_c')
```

### 对数刻度

```python
fig = px.scatter(df, x='x', y='y', log_x=True, log_y=True)
```

### 图像显示

```python
import plotly.express as px
fig = px.imshow(img_array)  # img_array from PIL, numpy, etc.
```

## 金融图表

### 蜡烛图

```python
import plotly.graph_objects as go
fig = go.Figure(data=[go.Candlestick(
    x=df['date'],
    open=df['open'],
    high=df['high'],
    low=df['low'],
    close=df['close']
)])
```

### OHLC 图表

```python
fig = go.Figure(data=[go.Ohlc(
    x=df['date'],
    open=df['open'],
    high=df['high'],
    low=df['low'],
    close=df['close']
)])
```

### 瀑布图

```python
fig = go.Figure(go.Waterfall(
    x=categories,
    y=values,
    measure=['relative', 'relative', 'total', 'relative', 'total']
))
```

### 漏斗图

```python
fig = px.funnel(df, x='count', y='stage')

# Or with graph_objects
fig = go.Figure(go.Funnel(
    y=['Stage 1', 'Stage 2', 'Stage 3'],
    x=[100, 60, 40]
))
```

### 时间序列

```python
fig = px.line(df, x='date', y='price')

# With rangeslider
fig.update_xaxes(rangeslider_visible=True)

# With range selector buttons
fig.update_xaxes(
    rangeselector=dict(
        buttons=list([
            dict(count=1, label='1m', step='month', stepmode='backward'),
            dict(count=6, label='6m', step='month', stepmode='backward'),
            dict(count=1, label='YTD', step='year', stepmode='todate'),
            dict(count=1, label='1y', step='year', stepmode='backward'),
            dict(step='all')
        ])
    )
)
```

## 地图和地理

### 散点图

```python
# Geographic projection
fig = px.scatter_geo(df, lat='lat', lon='lon', color='value', size='size')

# Mapbox (requires token for some styles)
fig = px.scatter_mapbox(
    df, lat='lat', lon='lon',
    color='value',
    zoom=10,
    mapbox_style='open-street-map'  # or 'carto-positron', 'carto-darkmatter'
)
```

### 等值线地图

```python
# Country-level
fig = px.choropleth(
    df,
    locations='iso_alpha',
    color='value',
    hover_name='country',
    color_continuous_scale='Viridis'
)

# US States
fig = px.choropleth(
    df,
    locations='state_code',
    locationmode='USA-states',
    color='value',
    scope='usa'
)
```

### 密度图

```python
fig = px.density_mapbox(
    df, lat='lat', lon='lon', z='value',
    radius=10,
    zoom=10,
    mapbox_style='open-street-map'
)
```

## 3D 图表

### 3D 散点图

```python
fig = px.scatter_3d(df, x='x', y='y', z='z', color='category', size='size')
```

### 3D 线

```python
fig = px.line_3d(df, x='x', y='y', z='z', color='group')
```

### 3D 表面

```python
import plotly.graph_objects as go
fig = go.Figure(data=[go.Surface(z=z_matrix, x=x_array, y=y_array)])

fig.update_layout(scene=dict(
    xaxis_title='X',
    yaxis_title='Y',
    zaxis_title='Z'
))
```

### 3D 网格

```python
fig = go.Figure(data=[go.Mesh3d(
    x=x_coords,
    y=y_coords,
    z=z_coords,
    i=i_indices,
    j=j_indices,
    k=k_indices,
    intensity=intensity_values,
    colorscale='Viridis'
)]
```

### 3D 锥体（矢量场）

```python
fig = go.Figure(data=go.Cone(
    x=x, y=y, z=z,
    u=u, v=v, w=w,
    colorscale='Blues',
    sizemode='absolute',
    sizeref=0.5
))
```

## 分层图表

### 旭日

```python
fig = px.sunburst(
    df,
    path=['continent', 'country', 'city'],
    values='population',
    color='value'
)
```

### 树状图

```python
fig = px.treemap(
    df,
    path=['category', 'subcategory', 'item'],
    values='count',
    color='value',
    color_continuous_scale='RdBu'
)
```

### 桑基图

```python
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color='black', width=0.5),
        label=['A', 'B', 'C', 'D', 'E'],
        color='blue'
    ),
    link=dict(
        source=[0, 1, 0, 2, 3],
        target=[2, 3, 3, 4, 4],
        value=[8, 4, 2, 8, 4]
    )
)])
```

## 专业图表

### 平行坐标

```python
fig = px.parallel_coordinates(
    df,
    dimensions=['dim1', 'dim2', 'dim3', 'dim4'],
    color='target',
    color_continuous_scale='Viridis'
)
```

### 并行类别

```python
fig = px.parallel_categories(
    df,
    dimensions=['cat1', 'cat2', 'cat3'],
    color='value'
)
```

### 散点矩阵 (SPLOM)

```python
fig = px.scatter_matrix(
    df,
    dimensions=['col1', 'col2', 'col3', 'col4'],
    color='category'
)
```

### 指示器/仪表

```python
fig = go.Figure(go.Indicator(
    mode='gauge+number+delta',
    value=75,
    delta={'reference': 60},
    gauge={'axis': {'range': [None, 100]},
           'bar': {'color': 'darkblue'},
           'steps': [
               {'range': [0, 50], 'color': 'lightgray'},
               {'range': [50, 100], 'color': 'gray'}
           ],
           'threshold': {'line': {'color': 'red', 'width': 4},
                        'thickness': 0.75,
                        'value': 90}
    }
))
```

### 表

```python
fig = go.Figure(data=[go.Table(
    header=dict(values=['A', 'B', 'C']),
    cells=dict(values=[col_a, col_b, col_c])
)])
```

## 生物信息学

### 树状图

```python
from plotly.figure_factory import create_dendrogram
fig = create_dendrogram(data_matrix)
```

### 带注释的热图

```python
from plotly.figure_factory import create_annotated_heatmap
fig = create_annotated_heatmap(z_matrix, x=x_labels, y=y_labels)
```

### 火山图

```python
# Typically built with scatter plot
fig = px.scatter(
    df,
    x='log2_fold_change',
    y='neg_log10_pvalue',
    color='significant',
    hover_data=['gene_name']
)
fig.add_hline(y=-np.log10(0.05), line_dash='dash')
fig.add_vline(x=-1, line_dash='dash')
fig.add_vline(x=1, line_dash='dash')
```