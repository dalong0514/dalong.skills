<!-- 此文件由机器翻译自 visualization.md -->

# 绘图和可视化

GeoPandas 通过 matplotlib 集成提供绘图。

## 基本绘图

```python
# Simple plot
gdf.plot()

# Customize figure size
gdf.plot(figsize=(10, 10))

# Set colors
gdf.plot(color='blue', edgecolor='black')

# Control line width
gdf.plot(edgecolor='black', linewidth=0.5)
```

## 等值线地图

基于数据值的颜色特征：

<<<代码块_1>>>

### 分类方案

需要：`uv pip install mapclassify`

<<<代码块_2>>>

### 图例定制

<<<代码块_3>>>

## 处理缺失数据

<<<代码块_4>>>

## 多层地图

组合多个 GeoDataFrame：

<<<代码块_5>>>

## 样式选项

<<<代码块_6>>>

## 地图增强

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12, 8))
gdf.plot(ax=ax, column='population', legend=True)

# Add title
ax.set_title('Population by Region', fontsize=16)

# Add axis labels
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Remove axes
ax.set_axis_off()

# Add north arrow and scale bar (requires separate packages)
# See geopandas-plot or contextily for these features

plt.tight_layout()
plt.show()
```

## 互动地图

需要：`uv pip install folium`

```python
# Create interactive map
m = gdf.explore(column='population', cmap='YlOrRd', legend=True)
m.save('map.html')

# Customize base map
m = gdf.explore(tiles='OpenStreetMap', legend=True)
m = gdf.explore(tiles='CartoDB positron', legend=True)

# Add tooltip
m = gdf.explore(column='population', tooltip=['name', 'population'], legend=True)

# Style options
m = gdf.explore(color='red', style_kwds={'fillOpacity': 0.5, 'weight': 2})

# Multiple layers
m = gdf1.explore(color='blue', name='Layer 1')
gdf2.explore(m=m, color='red', name='Layer 2')
folium.LayerControl().add_to(m)
```

## 与其他绘图类型集成

GeoPandas 支持 pandas 绘图类型：

```python
# Histogram of attribute
gdf['population'].plot.hist(bins=20)

# Scatter plot
gdf.plot.scatter(x='income', y='population')

# Box plot
gdf.boxplot(column='population', by='region')
```

## 带有上下文的底图

需要：`uv pip install contextily`

```python
import contextily as ctx

# Reproject to Web Mercator for basemap compatibility
gdf_webmercator = gdf.to_crs(epsg=3857)

fig, ax = plt.subplots(figsize=(10, 10))
gdf_webmercator.plot(ax=ax, alpha=0.5, edgecolor='k')

# Add basemap
ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)
# Other sources: ctx.providers.CartoDB.Positron, ctx.providers.Stamen.Terrain

plt.show()
```

## 使用 CartoPy 进行制图投影

需要：`uv pip install cartopy`

```python
import cartopy.crs as ccrs

# Create map with specific projection
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Robinson()}, figsize=(15, 10))

gdf.plot(ax=ax, transform=ccrs.PlateCarree(), column='population', legend=True)

ax.coastlines()
ax.gridlines(draw_labels=True)

plt.show()
```

## 保存数字

```python
# Save to file
ax = gdf.plot()
fig = ax.get_figure()
fig.savefig('map.png', dpi=300, bbox_inches='tight')
fig.savefig('map.pdf')
fig.savefig('map.svg')
```