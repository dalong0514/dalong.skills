<!-- 此文件由机器翻译自 data-io.md -->

# 读写空间数据

## 读取文件

使用 `geopandas.read_file()` 导入矢量空间数据：

```python
import geopandas as gpd

# Read from file
gdf = gpd.read_file("data.shp")
gdf = gpd.read_file("data.geojson")
gdf = gpd.read_file("data.gpkg")

# Read from URL
gdf = gpd.read_file("https://example.com/data.geojson")

# Read from ZIP archive
gdf = gpd.read_file("data.zip")
```

### 性能：箭头加速

要使阅读速度加快 2-4 倍，请使用箭头：

<<<代码块_1>>>

需要 PyArrow：`uv pip install pyarrow`

### 读取期间过滤

预过滤数据以仅加载需要的数据：

<<<代码块_2>>>

## 写入文件

使用 `to_file()` 导出：

<<<代码块_3>>>

### 支持的格式

列出所有可用的驱动程序：

<<<代码块_4>>>

常见格式：Shapefile、GeoJSON、GeoPackage (GPKG)、KML、MapInfo 文件、CSV（带有 WKT 几何图形）

## 镶木地板和羽毛

保留空间信息并支持多个几何列的列格式：

<<<代码块_5>>>

优点：
- 比传统格式更快的 I/O
- 更好的压缩
- 保留多个几何列
- 模式版本控制支持

## PostGIS 数据库

### 从 PostGIS 读取

<<<代码块_6>>>

### 写入 PostGIS

```python
# Create or replace table
gdf.to_postgis("table_name", con=engine, if_exists='replace')

# Append to existing table
gdf.to_postgis("table_name", con=engine, if_exists='append')

# Fail if table exists
gdf.to_postgis("table_name", con=engine, if_exists='fail')
```

需要：`uv pip install psycopg2` 或 `uv pip install psycopg` 和 `uv pip install geoalchemy2`

## 类文件对象

从文件句柄或内存缓冲区中读取：

```python
# From file handle
with open('data.geojson', 'r') as f:
    gdf = gpd.read_file(f)

# From StringIO
from io import StringIO
geojson_string = '{"type": "FeatureCollection", ...}'
gdf = gpd.read_file(StringIO(geojson_string))
```

## 远程存储（fsspec）

从云存储访问数据：

```python
# S3
gdf = gpd.read_file("s3://bucket/data.gpkg")

# Azure Blob Storage
gdf = gpd.read_file("az://container/data.gpkg")

# HTTP/HTTPS
gdf = gpd.read_file("https://example.com/data.geojson")
```