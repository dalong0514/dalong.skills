<!-- 此文件由机器翻译自 spatial-analysis.md -->

# 空间分析

## 属性连接

使用标准 pandas 合并基于公共变量合并数据集：

```python
# Merge on common column
result = gdf.merge(df, on='common_column')

# Left join
result = gdf.merge(df, on='common_column', how='left')

# Important: Call merge on GeoDataFrame to preserve geometry
# This works: gdf.merge(df, ...)
# This doesn't: df.merge(gdf, ...) # Returns DataFrame, not GeoDataFrame
```

## 空间连接

根据空间关系组合数据集。

### 二元谓词连接 (sjoin)

基于几何谓词的连接：

<<<代码块_1>>>

`how` 参数确定保留哪些几何图形：
- **左**：保留左GeoDataFrame的索引和几何图形
- **右**：保留正确的 GeoDataFrame 的索引和几何图形
- **内部**：使用索引的交集，保留左侧几何图形

### 最近连接 (sjoin_nearest)

加入最近的特征：

<<<代码块_2>>>

## 叠加操作

组合来自两个 GeoDataFrame 的几何图形的集合论运算：

<<<代码块_3>>>

结果包括来自两个输入 GeoDataFrame 的属性。

## 溶解（聚合）

根据属性值聚合几何图形：

<<<代码块_4>>>

## 剪辑

将几何图形裁剪到另一个几何图形的边界：

<<<代码块_5>>>

## 追加

组合多个 GeoDataFrame：

<<<代码块_6>>>

## 空间索引

提高空间操作的性能：

```python
# GeoPandas uses spatial index automatically for most operations
# Access the spatial index directly
sindex = gdf.sindex

# Query geometries intersecting a bounding box
possible_matches_index = list(sindex.intersection((xmin, ymin, xmax, ymax)))
possible_matches = gdf.iloc[possible_matches_index]

# Query geometries intersecting a polygon
possible_matches_index = list(sindex.query(polygon_geometry))
possible_matches = gdf.iloc[possible_matches_index]
```

空间索引显着加快：
- 空间连接
- 叠加操作
- 带有几何谓词的查询

## 距离计算

```python
# Distance between geometries
distances = gdf1.geometry.distance(gdf2.geometry)

# Distance to single geometry
distances = gdf.geometry.distance(single_point)

# Minimum distance to any feature
min_dist = gdf.geometry.distance(point).min()
```

## 面积和长度计算

为了准确测量，请确保正确的 CRS：

```python
# Reproject to appropriate projected CRS for area/length calculations
gdf_projected = gdf.to_crs(epsg=3857)  # Or appropriate UTM zone

# Calculate area (in CRS units, typically square meters)
areas = gdf_projected.geometry.area

# Calculate length/perimeter (in CRS units)
lengths = gdf_projected.geometry.length
```