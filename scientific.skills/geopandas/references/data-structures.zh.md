<!-- 此文件由机器翻译自 data-structures.md -->

# GeoPandas 数据结构

## 地理系列

GeoSeries 是一个向量，其中每个条目都是与一个观察相对应的一组形状（类似于 pandas Series，但具有几何数据）。

```python
import geopandas as gpd
from shapely.geometry import Point, Polygon

# Create a GeoSeries from geometries
points = gpd.GeoSeries([Point(1, 1), Point(2, 2), Point(3, 3)])

# Access geometric properties
points.area
points.length
points.bounds
```

## 地理数据框

GeoDataFrame 是包含 GeoSeries 的表格数据结构（类似于 pandas DataFrame，但包含地理数据）。

<<<代码块_1>>>

## 关键属性

- **geometry**：活动几何列（可以有多个几何列）
- **crs**：坐标参考系
- **bounds**：所有几何图形的边界框
- **total_bounds**：总体边界框

## 设置活动几何体

当 GeoDataFrame 有多个几何列时：

<<<代码块_2>>>

## 索引和选择

对空间数据使用标准 pandas 索引：

<<<代码块_3>>>