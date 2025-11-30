<!-- 此文件由机器翻译自 crs-management.md -->

# 坐标参考系统 (CRS)

坐标参考系统定义坐标与地球上位置的关系。

## 了解 CRS

CRS 信息存储为 `pyproj.CRS` 对象：

```python
# Check CRS
print(gdf.crs)

# Check if CRS is set
if gdf.crs is None:
    print("No CRS defined")
```

## 设置与重新投影

### 设置 CRS

当坐标正确但缺少 CRS 元数据时，使用 `set_crs()`：

<<<代码块_1>>>

**警告**：仅在 CRS 元数据丢失时使用。这不会转换坐标。

### 重新投影

使用 `to_crs()` 在坐标系之间转换坐标：

<<<代码块_2>>>

## CRS 格式

GeoPandas 通过 `pyproj.CRS.from_user_input()` 接受多种格式：

<<<代码块_3>>>

**最佳实践**：使用 WKT2 或权限字符串 (EPSG) 保留完整的 CRS 信息。

## 常见 EPSG 代码

### 地理坐标系

<<<代码块_4>>>

### 投影坐标系

<<<代码块_5>>>

## CRS 运营要求

### 需要匹配 CRS 的操作

这些操作需要相同的 CRS：

<<<代码块_6>>>

### 预计 CRS 中最佳运营

面积和距离计算应使用投影 CRS：

```python
# Bad: area in degrees (meaningless)
areas_degrees = gdf.geometry.area  # If CRS is EPSG:4326

# Good: reproject to appropriate projected CRS first
gdf_projected = gdf.to_crs("EPSG:3857")
areas_meters = gdf_projected.geometry.area  # Square meters

# Better: use appropriate local UTM zone for accuracy
gdf_utm = gdf.to_crs("EPSG:32633")  # UTM Zone 33N
accurate_areas = gdf_utm.geometry.area
```

## 选择合适的 CRS

### 用于面积/距离计算

使用等面积投影：

```python
# Albers Equal Area Conic (North America)
gdf.to_crs("EPSG:5070")

# Lambert Azimuthal Equal Area
gdf.to_crs("EPSG:3035")  # Europe

# UTM zones (for local areas)
gdf.to_crs("EPSG:32633")  # Appropriate UTM zone
```

### 用于保持距离（导航）

使用等距投影：

```python
# Azimuthal Equidistant
gdf.to_crs("ESRI:54032")
```

### 用于保持形状（角度）

使用等角投影：

```python
# Web Mercator (conformal but distorts area)
gdf.to_crs("EPSG:3857")

# UTM zones (conformal for local areas)
gdf.to_crs("EPSG:32633")
```

### 对于网络地图

```python
# Web Mercator (standard for web maps)
gdf.to_crs("EPSG:3857")
```

## 估计 UTM 区域

```python
# Estimate appropriate UTM CRS from data
utm_crs = gdf.estimate_utm_crs()
gdf_utm = gdf.to_crs(utm_crs)
```

## 具有不同 CRS 的多个几何列

GeoPandas 0.8+ 支持每个几何列不同的 CRS：

```python
# Set CRS for specific geometry column
gdf = gdf.set_crs("EPSG:4326", allow_override=True)

# Active geometry determines operations
gdf = gdf.set_geometry('other_geom_column')

# Check CRS mismatch
try:
    result = gdf1.overlay(gdf2)
except ValueError as e:
    print("CRS mismatch:", e)
```

## CRS 信息

```python
# Get full CRS details
print(gdf.crs)

# Get EPSG code if available
print(gdf.crs.to_epsg())

# Get WKT representation
print(gdf.crs.to_wkt())

# Get PROJ string
print(gdf.crs.to_proj4())

# Check if CRS is geographic (lat/lon)
print(gdf.crs.is_geographic)

# Check if CRS is projected
print(gdf.crs.is_projected)
```

## 变换单个几何形状

```python
from pyproj import Transformer

# Create transformer
transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)

# Transform point
x_new, y_new = transformer.transform(x, y)
```