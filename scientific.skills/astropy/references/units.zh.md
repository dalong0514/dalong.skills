<!-- 此文件由机器翻译自 units.md -->

# 单位和数量 (astropy.units)

`astropy.units` 模块处理物理量的定义、转换以及执行算术运算。

## 创建数量

将数值乘以或除以内置单位以创建 Quantity 对象：

```python
from astropy import units as u
import numpy as np

# Scalar quantities
distance = 42.0 * u.meter
velocity = 100 * u.km / u.s

# Array quantities
distances = np.array([1., 2., 3.]) * u.m
wavelengths = [500, 600, 700] * u.nm
```

通过 `.value` 和 `.unit` 属性访问组件：
<<<代码块_1>>>

## 单位换算

使用 `.to()` 方法进行转换：

<<<代码块_2>>>

## 算术运算

数量支持标准算术和自动单位管理：

<<<代码块_3>>>

## 单位系统

主要单位制之间的转换：

<<<代码块_4>>>

## 等价物

特定于域的转换需要等效项：

<<<代码块_5>>>

## 对数单位

震级、分贝和 dex 的特殊单位：

<<<代码块_6>>>

## 常用单位

### 长度
`u.m, u.km, u.cm, u.mm, u.micron, u.angstrom, u.au, u.pc, u.kpc, u.Mpc, u.lyr`

### 时间
`u.s, u.min, u.hour, u.day, u.year, u.Myr, u.Gyr`

### 弥撒
`u.kg, u.g, u.M_sun, u.M_earth, u.M_jup`

### 温度
`u.K, u.deg_C`

### 角度
`u.deg, u.arcmin, u.arcsec, u.rad, u.hourangle, u.mas`

### 能源/电力
`u.J, u.erg, u.eV, u.keV, u.MeV, u.GeV, u.W, u.L_sun`

### 频率
`u.Hz, u.kHz, u.MHz, u.GHz`

### 通量
`u.Jy, u.mJy, u.erg / u.s / u.cm**2`

## 性能优化

预先计算数组操作的复合单元：

```python
# Slow (creates intermediate quantities)
result = array * u.m / u.s / u.kg / u.sr

# Fast (pre-computed composite unit)
UNIT_COMPOSITE = u.m / u.s / u.kg / u.sr
result = array * UNIT_COMPOSITE

# Fastest (avoid copying with <<)
result = array << UNIT_COMPOSITE  # 10000x faster
```

## 字符串格式化

使用标准 Python 语法格式化数量：

```python
velocity = 15.1 * u.meter / (32.0 * u.second)
f"{velocity:0.03f}"     # '0.472 m / s'
f"{velocity:.2e}"       # '4.72e-01 m / s'
f"{velocity.unit:FITS}" # 'm s-1'
```

## 定义自定义单位

```python
# Create new unit
bakers_fortnight = u.def_unit('bakers_fortnight', 13 * u.day)

# Enable in string parsing
u.add_enabled_units([bakers_fortnight])
```

## 常量

访问带单位的物理常数：

```python
from astropy.constants import c, G, M_sun, h, k_B

speed_of_light = c.to(u.km/u.s)
gravitational_constant = G.to(u.m**3 / u.kg / u.s**2)
```