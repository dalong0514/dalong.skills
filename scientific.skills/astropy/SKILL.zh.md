<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 天体视学
描述：用于天文学和天体物理学的综合 Python 库。在处理天文数据（包括天体坐标、物理单位、FITS 文件、宇宙学计算、时间系统、表格、世界坐标系 (WCS) 和天文数据分析）时，应使用此技能。当任务涉及坐标变换、单位转换、FITS 文件操作、宇宙距离计算、时间尺度转换或天文数据处理时使用。
---

# 天体

## 概述

Astropy 是天文学的核心 Python 包，为天文学研究和数据分析提供基本功能。使用 astropy 进行坐标变换、单位和数量计算、FITS 文件操作、宇宙学计算、精确时间处理、表格数据操作和天文图像处理。

## 何时使用此技能

当任务涉及以下内容时使用天文学：
- 天球坐标系之间的转换（ICRS、Gactic、FK5、AltAz 等）
- 使用物理单位和数量（将 Jy 转换为 mJy、将秒差距转换为 km 等）
- 读取、写入或操作 FITS 文件（图像或表格）
- 宇宙学计算（光度距离、回溯时间、哈勃参数）
- 不同时间尺度（UTC、TAI、TT、TDB）和格式（JD、MJD、ISO）的精确时间处理
- 表操作（读取目录、交叉匹配、过滤、连接）
- 像素和世界坐标之间的 WCS 转换
- 天文常数和计算

## 快速入门

```python
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import Planck18

# Units and quantities
distance = 100 * u.pc
distance_km = distance.to(u.km)

# Coordinates
coord = SkyCoord(ra=10.5*u.degree, dec=41.2*u.degree, frame='icrs')
coord_galactic = coord.galactic

# Time
t = Time('2023-01-15 12:30:00')
jd = t.jd  # Julian Date

# FITS files
data = fits.getdata('image.fits')
header = fits.getheader('image.fits')

# Tables
table = Table.read('catalog.fits')

# Cosmology
d_L = Planck18.luminosity_distance(z=1.0)
```

## 核心能力

### 1. 单位和数量 (`astropy.units`)

处理带有单位的物理量，执行单位转换，并确保计算中的量纲一致性。

**关键操作：**
- 通过将值与单位相乘来创建数量
- 使用 `.to()` 方法在单位之间进行转换
- 通过自动单位处理执行算术
- 使用等效项进行特定领域的转换（光谱、多普勒、视差）
- 使用对数单位（幅度、分贝）

**请参阅：** `references/units.md` 了解全面的文档、单位系统、等效项、性能优化和单位算术。

### 2. 坐标系 (`astropy.coordinates`)

表示天体位置并在不同坐标系之间进行变换。

**关键操作：**
- 在任何帧（ICRS、Gactic、FK5、AltAz 等）中使用 `SkyCoord` 创建坐标
- 坐标系之间的变换
- 计算角距和位置角
- 将坐标与目录相匹配
- 包括 3D 坐标操作的距离
- 处理自行和径向速度
- 从在线数据库查询命名对象

**请参阅：** `references/coordinates.md` 了解详细的坐标系描述、转换、依赖于观察者的坐标系 (AltAz)、目录匹配和性能提示。

### 3. 宇宙学计算 (`astropy.cosmology`)

使用标准宇宙学模型执行宇宙学计算。

**关键操作：**
- 使用内置宇宙学（Planck18、WMAP9 等）
- 创建自定义宇宙模型
- 计算距离（光度、同移、角直径）
- 计算年龄和回顾时间
- 确定任意红移处的哈勃参数
- 计算密度参数和体积
- 执行逆计算（找到给定距离的 z）

**请参阅：** `references/cosmology.md` 了解可用模型、距离计算、时间计算、密度参数和中微子效应。

### 4. FITS 文件处理 (`astropy.io.fits`)

读取、写入和操作 FITS（灵活图像传输系统）文件。

**关键操作：**
- 使用上下文管理器打开 FITS 文件
- 按索引或名称访问 HDU（标头数据单元）
- 读取和修改标题（关键词、评论、历史记录）
- 处理图像数据（NumPy 数组）
- 处理表数据（二进制和 ASCII 表）
- 创建新的 FITS 文件（单个或多个扩展名）
- 对大文件使用内存映射
- 访问远程 FITS 文件（S3、HTTP）

**请参阅：** `references/fits.md` 了解全面的文件操作、标头操作、图像和表处理、多扩展文件以及性能注意事项。

### 5. 表操作 (`astropy.table`)

处理表格数据并支持单位、元数据和各种文件格式。

**关键操作：**
- 从数组、列表或字典创建表
- 读/写多种格式的表（FITS、CSV、HDF5、VOTable）
- 访问和修改列和行
- 对表进行排序、过滤和索引
- 执行数据库式操作（连接、分组、聚合）
- 堆叠和连接表格
- 使用单位感知列（QTable）
- 通过屏蔽处理缺失数据

**请参阅：** `references/tables.md` 了解表创建、I/O 操作、数据操作、排序、过滤、连接、分组和性能提示。

### 6. 时间处理 (`astropy.time`)

精确的时间表示以及时间尺度和格式之间的转换。

**关键操作：**
- 创建各种格式的时间对象（ISO、JD、MJD、Unix 等）
- 时间尺度之间的转换（UTC、TAI、TT、TDB 等）
- 使用 TimeDelta 执行时间算术
- 计算观测者的恒星时
- 计算光传播时间修正（重心、日心）
- 有效地处理时间数组
- 处理屏蔽（缺失）时间

**请参阅：** `references/time.md` 了解时间格式、时间尺度、转换、算术、观察功能和精度处理。

### 7. 世界坐标系 (`astropy.wcs`)

图像中的像素坐标和世界坐标之间的转换。

**关键操作：**
- 从 FITS 标头读取 WCS
- 将像素坐标转换为世界坐标（反之亦然）
- 计算图像足迹
- 访问 WCS 参数（参考像素、投影、比例）
- 创建自定义 WCS 对象

**请参阅：** `references/wcs_and_other_modules.md` 了解 WCS 操作和转换。

## 附加功能

`references/wcs_and_other_modules.md` 文件还涵盖：

### NDData 和 CCDData
包含元数据、不确定性、屏蔽和 WCS 信息的 n 维数据集的容器。

### 建模
用于创建数学模型并将其拟合到天文数据的框架。

### 可视化
通过适当的拉伸和缩放来显示天文图像的工具。

### 常量
具有适当单位的物理和天文常数（光速、太阳质量、普朗克常数等）。

### 卷积
用于平滑和过滤的图像处理内核。

### 统计
强大的统计功能，包括西格玛剪裁和异常值拒绝。

## 安装

<<<代码块_1>>>

## 常见工作流程

### 在系统之间转换坐标

<<<代码块_2>>>

### 读取和分析 FITS 文件

<<<代码块_3>>>

### 宇宙距离计算

<<<代码块_4>>>

### 交叉匹配目录

<<<代码块_5>>>

## 最佳实践

1. **始终使用单位**：将单位附加到数量以避免错误并确保尺寸一致性
2. **对 FITS 文件使用上下文管理器**：确保正确的文件关闭
3. **优先使用数组而不是循环**：将多个坐标/时间作为数组处理以获得更好的性能
4. **检查坐标系**：在变换之前验证坐标系
5. **使用适当的宇宙学**：为您的分析选择正确的宇宙学模型
6. **处理缺失数据**：对缺失值的表使用屏蔽列
7. **指定时间尺度**：明确时间尺度（UTC、TT、TDB）以实现精确计时
8. **将 QTable 用于单位感知表**：当表列具有单位时
9. **检查 WCS 有效性**：在使用转换之前验证 WCS
10. **缓存常用值**：可以缓存昂贵的计算（例如宇宙距离）

## 文档和资源

- Astropy 官方文档：https://docs.astropy.org/en/stable/
- 教程：https://learn.astropy.org/
- GitHub：https://github.com/astropy/astropy

## 参考文件

有关特定模块的详细信息：
- `references/units.md` - 单位、数量、换算和等价
- `references/coordinates.md` - 坐标系、变换和目录匹配
- `references/cosmology.md` - 宇宙学模型和计算
- `references/fits.md` - 适合文件操作和操作
- `references/tables.md` - 表创建、I/O 和操作
- `references/time.md` - 时间格式、比例和计算
- `references/wcs_and_other_modules.md` - WCS、NDData、建模、可视化、常量和实用程序