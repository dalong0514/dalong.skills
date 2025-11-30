<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：地理熊猫
描述：用于处理地理空间矢量数据（包括 shapefile、GeoJSON 和 GeoPackage 文件）的 Python 库。在处理地理数据进行空间分析、几何运算、坐标变换、空间连接、叠加操作、分区统计图或涉及读取/写入/分析矢量地理数据的任何任务时使用。支持 PostGIS 数据库、交互式地图以及与 matplotlib/folium/cartopy 的集成。用于缓冲区分析、数据集之间的空间连接、溶解边界、裁剪数据、计算面积/距离、重新投影坐标系、创建地图或在空间文件格式之间进行转换等任务。
---

# GeoPandas

GeoPandas 扩展了 pandas 以支持几何类型的空间操作。它结合了 pandas 和 shapely 的地理空间数据分析功能。

## 安装

```bash
uv pip install geopandas
```

### 可选依赖项

<<<代码块_1>>>

## 快速入门

<<<代码块_2>>>

## 核心概念

### 数据结构

- **GeoSeries**：具有空间操作的几何向量
- **GeoDataFrame**：具有几何列的表格数据结构

详细信息请参见[data-structs.md](references/data-structurals.md)。

### 读取和写入数据

GeoPandas 读取/写入多种格式：Shapefile、GeoJSON、GeoPackage、PostGIS、Parquet。

<<<代码块_3>>>

有关全面的 I/O 操作，请参阅 [data-io.md](references/data-io.md)。

### 坐标参考系

始终检查和管理 CRS 以实现准确的空间操作：

<<<代码块_4>>>

CRS操作请参见[crs-management.md](references/crs-management.md)。

## 常用操作

### 几何运算

缓冲、简化、质心、凸包、仿射变换：

<<<代码块_5>>>

有关所有操作，请参阅 [geometric-operations.md](references/geometric-operations.md)。

### 空间分析

空间连接、叠加操作、溶解：

<<<代码块_6>>>

分析操作参见[spatial-analysis.md](references/spatial-analysis.md)。

### 可视化

创建静态和交互式地图：

```python
# Choropleth map
gdf.plot(column='population', cmap='YlOrRd', legend=True)

# Interactive map
gdf.explore(column='population', legend=True).save('map.html')

# Multi-layer map
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
gdf1.plot(ax=ax, color='blue')
gdf2.plot(ax=ax, color='red')
```

有关映射技术，请参阅 [visualization.md](references/visualization.md)。

## 详细文档

- **[数据结构](references/data-structs.md)** - GeoSeries 和 GeoDataFrame 基础知识
- **[数据 I/O](references/data-io.md)** - 读/写文件、PostGIS、Parquet
- **[几何操作](references/geometric-operations.md)** - 缓冲、简化、仿射变换
- **[空间分析](references/spatial-analysis.md)** - 连接、叠加、溶解、裁剪
- **[可视化](references/visualization.md)** - 绘图、分区统计图、交互式地图
- **[CRS Management](references/crs-management.md)** - 坐标参考系统和投影

## 常见工作流程

### 加载、转换、分析、导出

```python
# 1. Load data
gdf = gpd.read_file("data.shp")

# 2. Check and transform CRS
print(gdf.crs)
gdf = gdf.to_crs("EPSG:3857")

# 3. Perform analysis
gdf['area'] = gdf.geometry.area
buffered = gdf.copy()
buffered['geometry'] = gdf.geometry.buffer(100)

# 4. Export results
gdf.to_file("results.gpkg", layer='original')
buffered.to_file("results.gpkg", layer='buffered')
```

### 空间连接和聚合

```python
# Join points to polygons
points_in_polygons = gpd.sjoin(points_gdf, polygons_gdf, predicate='within')

# Aggregate by polygon
aggregated = points_in_polygons.groupby('index_right').agg({
    'value': 'sum',
    'count': 'size'
})

# Merge back to polygons
result = polygons_gdf.merge(aggregated, left_index=True, right_index=True)
```

### 多源数据集成

```python
# Read from different sources
roads = gpd.read_file("roads.shp")
buildings = gpd.read_file("buildings.geojson")
parcels = gpd.read_postgis("SELECT * FROM parcels", con=engine, geom_col='geom')

# Ensure matching CRS
buildings = buildings.to_crs(roads.crs)
parcels = parcels.to_crs(roads.crs)

# Perform spatial operations
buildings_near_roads = buildings[buildings.geometry.distance(roads.union_all()) < 50]
```

## 性能提示

1. **使用空间索引**：GeoPandas 会为大多数操作自动创建空间索引
2. **读取期间过滤**：使用`bbox`、`mask`或`where`参数仅加载需要的数据
3. **使用 Arrow 进行 I/O**：添加 `use_arrow=True` 使读/写速度提高 2-4 倍
4. **简化几何图形**：当精度不重要时，使用 `.simplify()` 来降低复杂性
5. **批量操作**：矢量化操作比迭代行快得多
6. **使用适当的 CRS**：面积/距离的投影 CRS，地理可视化

## 最佳实践

1. **在空间操作之前始终检查 CRS**
2. **使用投影 CRS** 进行面积和距离计算
3. **在空间连接或覆盖之前匹配 CRS**
4. **在操作之前使用 `.is_valid` 验证几何图形**
5. **修改几何列时使用`.copy()`**以避免副作用
6. **在简化分析时保留拓扑**
7. **使用 GeoPackage** 格式进行现代工作流程（比 Shapefile 更好）
8. **在 sjoin_nearest 中设置 max_distance** 以获得更好的性能