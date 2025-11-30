<!-- 此文件由机器翻译自 geometric-operations.md -->

# 几何运算

GeoPandas 通过 Shapely 集成提供广泛的几何操作。

## 建设性操作

从现有几何图形创建新的几何图形：

### 缓冲区

创建代表一定距离内所有点的几何图形：

```python
# Buffer by fixed distance
buffered = gdf.geometry.buffer(10)

# Negative buffer (erosion)
eroded = gdf.geometry.buffer(-5)

# Buffer with resolution parameter
smooth_buffer = gdf.geometry.buffer(10, resolution=16)
```

### 边界

获取低维边界：

<<<代码块_1>>>

### 质心

获取每个几何体的中心点：

<<<代码块_2>>>

### 凸包

包含所有点的最小凸多边形：

<<<代码块_3>>>

### 凹形船体

包含所有点的最小凹多边形：

<<<代码块_4>>>

### 信封

最小的轴对齐矩形：

<<<代码块_5>>>

### 简化

降低几何复杂度：

<<<代码块_6>>>

### 细分

向线段添加顶点：

```python
# Add vertices with maximum segment length
segmented = gdf.geometry.segmentize(max_segment_length=5)
```

### 联合所有

将所有几何图形合并为单个几何图形：

```python
# Union all features
unified = gdf.geometry.union_all()
```

## 仿射变换

坐标的数学变换：

### 旋转

```python
# Rotate around origin (0, 0) by angle in degrees
rotated = gdf.geometry.rotate(angle=45, origin='center')

# Rotate around custom point
rotated = gdf.geometry.rotate(angle=45, origin=(100, 100))
```

### 规模

```python
# Scale uniformly
scaled = gdf.geometry.scale(xfact=2.0, yfact=2.0)

# Scale with origin
scaled = gdf.geometry.scale(xfact=2.0, yfact=2.0, origin='center')
```

### 翻译

```python
# Shift coordinates
translated = gdf.geometry.translate(xoff=100, yoff=50)
```

### 倾斜

```python
# Shear transformation
skewed = gdf.geometry.skew(xs=15, ys=0, origin='center')
```

### 自定义仿射变换

```python
from shapely import affinity

# Apply 6-parameter affine transformation matrix
# [a, b, d, e, xoff, yoff]
transformed = gdf.geometry.affine_transform([1, 0, 0, 1, 100, 50])
```

## 几何属性

访问几何属性（返回 pandas 系列）：

```python
# Area
areas = gdf.geometry.area

# Length/perimeter
lengths = gdf.geometry.length

# Bounding box coordinates
bounds = gdf.geometry.bounds  # Returns DataFrame with minx, miny, maxx, maxy

# Total bounds for entire GeoSeries
total_bounds = gdf.geometry.total_bounds  # Returns array [minx, miny, maxx, maxy]

# Check geometry types
geom_types = gdf.geometry.geom_type

# Check if valid
is_valid = gdf.geometry.is_valid

# Check if empty
is_empty = gdf.geometry.is_empty
```

## 几何关系

二元谓词测试关系：

```python
# Within
gdf1.geometry.within(gdf2.geometry)

# Contains
gdf1.geometry.contains(gdf2.geometry)

# Intersects
gdf1.geometry.intersects(gdf2.geometry)

# Touches
gdf1.geometry.touches(gdf2.geometry)

# Crosses
gdf1.geometry.crosses(gdf2.geometry)

# Overlaps
gdf1.geometry.overlaps(gdf2.geometry)

# Covers
gdf1.geometry.covers(gdf2.geometry)

# Covered by
gdf1.geometry.covered_by(gdf2.geometry)
```

## 点提取

从几何图形中提取特定点：

```python
# Representative point (guaranteed to be within geometry)
rep_points = gdf.geometry.representative_point()

# Interpolate point along line at distance
points = line_gdf.geometry.interpolate(distance=10)

# Interpolate point at normalized distance (0 to 1)
midpoints = line_gdf.geometry.interpolate(distance=0.5, normalized=True)
```

## Delaunay 三角测量

```python
# Create triangulation
triangles = gdf.geometry.delaunay_triangles()
```