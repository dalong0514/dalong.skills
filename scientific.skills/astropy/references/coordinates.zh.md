<!-- 此文件由机器翻译自 coordinates.md -->

# 天文坐标 (astropy.coordinates)

`astropy.coordinates` 包提供了用于表示天体坐标并在不同坐标系之间进行转换的工具。

## 使用 SkyCoord 创建坐标

高级 `SkyCoord` 类是推荐的接口：

```python
from astropy import units as u
from astropy.coordinates import SkyCoord

# Decimal degrees
c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')

# Sexagesimal strings
c = SkyCoord(ra='00h42m30s', dec='+41d12m00s', frame='icrs')

# Mixed formats
c = SkyCoord('00h42.5m +41d12m', unit=(u.hourangle, u.deg))

# Galactic coordinates
c = SkyCoord(l=120.5*u.degree, b=-23.4*u.degree, frame='galactic')
```

## 数组坐标

使用数组有效处理多个坐标：

<<<代码块_1>>>

## 访问组件

<<<代码块_2>>>

## 字符串格式化

<<<代码块_3>>>

## 坐标变换

参考系之间的变换：

<<<代码块_4>>>

## 常用坐标系

### 天体框架
- **ICRS**：国际天体参考系统（默认，最常见）
- **FK5**：第五基本目录（默认为Equinox J2000.0）
- **FK4**：第四个基本目录（较旧，需要春分规范）
- **GCRS**：地心天体参考系统
- **CIRS**：天体中间参考系统

### 银河框架
- **银河**：IAU 1958 银河坐标
- **超级银河**：De Vaucouleurs 超级银河坐标
- **银河中心**：基于银河中心的 3D 坐标

### 水平框架
- **AltAz**：高度-方位角（取决于观察者）
- **HADec**：时角偏角

### 黄道框架
- **地心平均黄道**：地心平均黄道
- **重心平均黄道**：重心平均黄道
- **日心平均黄道**：日心平均黄道

## 观察者相关的变换

对于高度-方位角坐标，指定观测时间和位置：

<<<代码块_5>>>

## 处理距离

添加 3D 坐标的距离信息：

<<<代码块_6>>>

## 角距

计算天空分离：

```python
c1 = SkyCoord(ra=10*u.degree, dec=9*u.degree, frame='icrs')
c2 = SkyCoord(ra=11*u.degree, dec=10*u.degree, frame='fk5')

# Angular separation (handles frame conversion automatically)
sep = c1.separation(c2)
print(f"Separation: {sep.arcsec} arcsec")

# Position angle
pa = c1.position_angle(c2)
```

## 目录匹配

将坐标与目录源匹配：

```python
# Single target matching
catalog = SkyCoord(ra=ra_array*u.degree, dec=dec_array*u.degree)
target = SkyCoord(ra=10.5*u.degree, dec=41.2*u.degree)

# Find closest match
idx, sep2d, dist3d = target.match_to_catalog_sky(catalog)
matched_coord = catalog[idx]

# Match with maximum separation constraint
matches = target.separation(catalog) < 1*u.arcsec
```

## 命名对象

从在线目录中检索坐标：

```python
# Query by name (requires internet)
m31 = SkyCoord.from_name("M31")
crab = SkyCoord.from_name("Crab Nebula")
psr = SkyCoord.from_name("PSR J1012+5307")
```

## 地球位置

定义观察者位置：

```python
# By coordinates
location = EarthLocation(lat=40*u.deg, lon=-120*u.deg, height=1000*u.m)

# By named observatory
keck = EarthLocation.of_site('Keck Observatory')
vlt = EarthLocation.of_site('Paranal Observatory')

# By address (requires internet)
location = EarthLocation.of_address('1002 Holy Grail Court, St. Louis, MO')

# List available observatories
EarthLocation.get_site_names()
```

## 速度信息

包括自行和径向速度：

```python
# Proper motion
c = SkyCoord(ra=10*u.degree, dec=41*u.degree,
             pm_ra_cosdec=15*u.mas/u.yr,
             pm_dec=5*u.mas/u.yr,
             distance=150*u.pc)

# Radial velocity
c = SkyCoord(ra=10*u.degree, dec=41*u.degree,
             radial_velocity=20*u.km/u.s)

# Both
c = SkyCoord(ra=10*u.degree, dec=41*u.degree, distance=150*u.pc,
             pm_ra_cosdec=15*u.mas/u.yr, pm_dec=5*u.mas/u.yr,
             radial_velocity=20*u.km/u.s)
```

## 表示类型

在坐标表示之间切换：

```python
# Cartesian representation
c = SkyCoord(x=1*u.kpc, y=2*u.kpc, z=3*u.kpc,
             representation_type='cartesian', frame='icrs')

# Change representation
c.representation_type = 'cylindrical'
c.rho  # Cylindrical radius
c.phi  # Azimuthal angle
c.z    # Height

# Spherical (default for most frames)
c.representation_type = 'spherical'
```

## 性能提示

1. **使用数组，而不是循环**：将多个坐标作为单个数组处理
2. **预计算帧**：重用帧对象进行多个转换
3. **使用广播**：多次高效地变换多个位置
4. **启用插值**：对于密集时间采样，请使用 ErfaAstromInterpolator

```python
# Fast approach
coords = SkyCoord(ra=ra_array*u.degree, dec=dec_array*u.degree)
coords_transformed = coords.transform_to('galactic')

# Slow approach (avoid)
for ra, dec in zip(ra_array, dec_array):
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    c_transformed = c.transform_to('galactic')
```