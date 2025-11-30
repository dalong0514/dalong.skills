<!-- 此文件由机器翻译自 core_classes.md -->

# Pymatgen 核心类参考

本参考文档记录了 `pymatgen.core` 中构成材料分析基础的基本类。

## 架构原则

Pymatgen 遵循面向对象的设计，其中元素、站点和结构都表示为对象。该框架强调晶体表示的周期性边界条件，同时保持分子系统的灵活性。

**单位约定**： pymatgen 中的所有单位通常假定为原子单位：
- 长度：埃 (Å)
- 能量：电子伏特 (eV)
- 角度：度

## 元素和元素周期表

### 元素
代表具有综合性质的元素周期表元素。

**创作方法：**
```python
from pymatgen.core import Element

# Create from symbol
si = Element("Si")
# Create from atomic number
si = Element.from_Z(14)
# Create from name
si = Element.from_name("silicon")
```

**关键属性：**
- `atomic_mass`：amu 中的原子质量
- `atomic_radius`：以埃为单位的原子半径
- `electronegativity`：鲍林电负性
- `ionization_energy`：第一电离能（eV）
- `common_oxidation_states`：常见氧化态列表
- `is_metal`、`is_halogen`、`is_noble_gas` 等：布尔属性
- `X`：元素符号作为字符串

### 物种
扩展带电离子和特定氧化态的元素。

<<<代码块_1>>>

### 虚拟物种
特殊结构表示的占位符原子（例如空位）。

<<<代码块_2>>>

## 成分

代表化学式和成分，可进行化学分析和操作。

### 创造
<<<代码块_3>>>

### 关键方法
- `get_reduced_formula_and_factor()`：返回简化公式和乘法因子
- `oxi_state_guesses()`：尝试确定氧化态
- `replace(replacements_dict)`：替换元素
- `add_charges_from_oxi_state_guesses()`：推断并添加氧化态
- `is_element`：检查组合是否是单个元素

### 关键属性
- `weight`：分子量
- `reduced_formula`：简化化学式
- `hill_formula`：希尔表示法的公式（C、H，然后按字母顺序）
- `num_atoms`：原子总数
- `chemical_system`：按字母顺序排序的元素（例如“Fe-O”）
- `element_composition`：要金额的元素字典

## 格子

定义晶体结构的晶胞几何形状。

### 创造
<<<代码块_4>>>

### 关键方法
- `get_niggli_reduced_lattice()`：返回 Niggli 简化格子
- `get_distance_and_image(frac_coords1, frac_coords2)`：具有周期性边界条件的分数坐标之间的距离
- `get_all_distances(frac_coords1, frac_coords2)`：包括周期性图像的距离

### 关键属性
- `volume`：晶胞体积 (Å³)
- `abc`：格子参数（a，b，c）作为元组
- `angles`：晶格角度（alpha、beta、gamma）作为元组
- `matrix`：3x3 晶格向量矩阵
- `reciprocal_lattice`：倒数晶格对象
- `is_orthogonal`：晶格向量是否正交

## 站点

### 网站
表示非周期系统中的原子位置。

<<<代码块_5>>>

### 定期网站
用分数坐标表示周期晶格中的原子位置。

<<<代码块_6>>>

**关键方法：**
- `distance(other_site)`：到另一个站点的距离
- `is_periodic_image(other_site)`：检查站点是否是周期性图像

**关键属性：**
- `species`：该地点的物种或元素
- `coords`：笛卡尔坐标
- `frac_coords`：分数坐标（对于PeriodicSite）
- `x`、`y`、`z`：单个笛卡尔坐标

## 结构

将晶体结构表示为周期性位点的集合。 `Structure` 是可变的，而 `IStructure` 是不可变的。

### 创造
```python
from pymatgen.core import Structure, Lattice

# From scratch
coords = [[0, 0, 0], [0.75, 0.5, 0.75]]
lattice = Lattice.from_parameters(a=3.84, b=3.84, c=3.84,
                                  alpha=120, beta=90, gamma=60)
struct = Structure(lattice, ["Si", "Si"], coords)

# From file (automatic format detection)
struct = Structure.from_file("POSCAR")
struct = Structure.from_file("structure.cif")

# From spacegroup
struct = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.5),
                                   ["Si"], [[0, 0, 0]])
```

### 文件输入/输出
```python
# Write to file (format inferred from extension)
struct.to(filename="output.cif")
struct.to(filename="POSCAR")
struct.to(filename="structure.xyz")

# Get string representation
cif_string = struct.to(fmt="cif")
poscar_string = struct.to(fmt="poscar")
```

### 关键方法

**结构修改：**
- `append(species, coords)`：添加站点
- `insert(i, species, coords)`：在索引处插入站点
- `remove_sites(indices)`：按索引删除站点
- `replace(i, species)`：替换索引处的物种
- `apply_strain(strain)`：对结构施加应变
- `perturb(distance)`：随机扰动原子位置
- `make_supercell(scaling_matrix)`：创建超级单元
- `get_primitive_structure()`：获取原始单元格

**分析：**
- `get_distance(i, j)`：站点 i 和 j 之间的距离
- `get_neighbors(site, r)`：获取半径 r 内的邻居
- `get_all_neighbors(r)`：获取所有站点的所有邻居
- `get_space_group_info()`：获取空间组信息
- `matches(other_struct)`：检查结构是否匹配

**插值：**
- `interpolate(end_structure, nimages)`：在结构之间插值
### 关键属性
- `lattice`：晶格对象
- `species`：每个地点的物种列表
- `sites`：PeriodicSite 对象列表
- `num_sites`：站点数量
- `volume`：结构体的体积
- `density`：密度，单位为 g/cm³
- `composition`：组合对象
- `formula`：化学式
- `distance_matrix`：成对距离矩阵

## 分子

表示原子的非周期性集合。 `Molecule` 是可变的，而 `IMolecule` 是不可变的。

### 创造
```python
from pymatgen.core import Molecule

# From scratch
coords = [[0.00, 0.00, 0.00],
          [0.00, 0.00, 1.08]]
mol = Molecule(["C", "O"], coords)

# From file
mol = Molecule.from_file("molecule.xyz")
mol = Molecule.from_file("molecule.mol")
```

### 关键方法
- `get_covalent_bonds()`：根据共价半径返回键
- `get_neighbors(site, r)`：获取半径内的邻居
- `get_zmatrix()`：获取 Z 矩阵表示
- `get_distance(i, j)`：站点之间的距离
- `get_centered_molecule()`：原点的中心分子

### 关键属性
- `species`：物种列表
- `sites`：站点对象列表
- `num_sites`：原子数量
- `charge`：分子的总电荷
- `spin_multiplicity`：自旋重数
- `center_of_mass`：质心坐标

## 序列化

所有核心对象都实现 `as_dict()` 和 `from_dict()` 方法以实现强大的 JSON/YAML 持久性。

```python
# Serialize to dictionary
struct_dict = struct.as_dict()

# Write to JSON
import json
with open("structure.json", "w") as f:
    json.dump(struct_dict, f)

# Read from JSON
with open("structure.json", "r") as f:
    struct_dict = json.load(f)
    struct = Structure.from_dict(struct_dict)
```

这种方法解决了 Python pickling 的限制并保持了 pymatgen 版本之间的兼容性。

## 其他核心课程

### 共价键
代表分子中的键。

**关键属性：**
- `length`：键长
- `get_bond_order()`：返回键顺序（单、双、三）

### 离子
代表具有氧化态的带电离子种类。

```python
from pymatgen.core import Ion

# Create Fe2+ ion
fe2_ion = Ion.from_formula("Fe2+")
```

### 接口
代表异质结分析的基底-薄膜组合。

### 晶界
代表晶体晶界。

### 频谱
使用归一化和处理方法表示光谱数据。

**关键方法：**
- `normalize(mode="max")`：标准化频谱
- `smear(sigma)`：应用高斯涂抹

## 最佳实践

1. **不可变性**：当不应修改结构时，使用不可变版本（`IStructure`、`IMolecule`）
2. **序列化**：对于长期存储，优先使用 `as_dict()`/`from_dict()` 而不是 pickle
3. **单位**：始终以原子单位（Å、eV）工作 - 可在 `pymatgen.core.units` 中进行转换
4. **文件I/O**：使用`from_file()`进行自动格式检测
5. **坐标**：注意方法是否期望笛卡尔坐标或分数坐标