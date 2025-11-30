<!-- 此文件由机器翻译自 io_formats.md -->

# Pymatgen I/O 和文件格式参考

该参考文档记录了 pymatgen 广泛的输入/输出功能，用于读取和写入 100 多种文件格式的结构和计算数据。

## 一般 I/O 原理

Pymatgen 通过 `from_file()` 和 `to()` 方法为文件操作提供统一的接口，并具有基于文件扩展名的自动格式检测。

### 读取文件

```python
from pymatgen.core import Structure, Molecule

# Automatic format detection
struct = Structure.from_file("POSCAR")
struct = Structure.from_file("structure.cif")
mol = Molecule.from_file("molecule.xyz")

# Explicit format specification
struct = Structure.from_file("file.txt", fmt="cif")
```

### 写入文件

<<<代码块_1>>>

## 结构文件格式

### CIF（晶体信息文件）
晶体学数据的标准格式。

<<<代码块_2>>>

**主要特点：**
- 支持对称信息
- 可以包含多个结构
- 保留空间群和对称操作
- 处理部分占用情况

### 波斯卡/CONTCAR (VASP)
VASP的结构格式。

<<<代码块_3>>>

**主要特点：**
- 支持选择性动态
- 可以包括速度（XDATCAR 格式）
- 保留晶格和坐标精度

### XYZ
简单的分子坐标格式。

<<<代码块_4>>>

### PDB（蛋白质数据库）
生物分子的通用格式。

<<<代码块_5>>>

### JSON/YAML
通过字典进行序列化。

<<<代码块_6>>>

## 电子结构代码I/O

### VASP

pymatgen 中最全面的集成。

#### 输入文件

```python
from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints, VaspInput

# INCAR (calculation parameters)
incar = Incar.from_file("INCAR")
incar = Incar({"ENCUT": 520, "ISMEAR": 0, "SIGMA": 0.05})
incar.write_file("INCAR")

# KPOINTS (k-point mesh)
from pymatgen.io.vasp.inputs import Kpoints
kpoints = Kpoints.automatic(20)  # 20x20x20 Gamma-centered mesh
kpoints = Kpoints.automatic_density(struct, 1000)  # By density
kpoints.write_file("KPOINTS")

# POTCAR (pseudopotentials)
potcar = Potcar(["Fe_pv", "O"])  # Specify functional variants

# Complete input set
vasp_input = VaspInput(incar, kpoints, poscar, potcar)
vasp_input.write_input("./vasp_calc")
```

#### 输出文件

```python
from pymatgen.io.vasp.outputs import Vasprun, Outcar, Oszicar, Eigenval

# vasprun.xml (comprehensive output)
vasprun = Vasprun("vasprun.xml")
final_structure = vasprun.final_structure
energy = vasprun.final_energy
band_structure = vasprun.get_band_structure()
dos = vasprun.complete_dos

# OUTCAR
outcar = Outcar("OUTCAR")
magnetization = outcar.total_mag
elastic_tensor = outcar.elastic_tensor

# OSZICAR (convergence information)
oszicar = Oszicar("OSZICAR")
```

#### 输入集

Pymatgen 为常见计算提供了预配置的输入集：

```python
from pymatgen.io.vasp.sets import (
    MPRelaxSet,      # Materials Project relaxation
    MPStaticSet,     # Static calculation
    MPNonSCFSet,     # Non-self-consistent (band structure)
    MPSOCSet,        # Spin-orbit coupling
    MPHSERelaxSet,   # HSE06 hybrid functional
)

# Create input set
relax = MPRelaxSet(struct)
relax.write_input("./relax_calc")

# Customize parameters
static = MPStaticSet(struct, user_incar_settings={"ENCUT": 600})
static.write_input("./static_calc")
```

### 高斯

量子化学包集成。

```python
from pymatgen.io.gaussian import GaussianInput, GaussianOutput

# Input
gin = GaussianInput(
    mol,
    charge=0,
    spin_multiplicity=1,
    functional="B3LYP",
    basis_set="6-31G(d)",
    route_parameters={"Opt": None, "Freq": None}
)
gin.write_file("input.gjf")

# Output
gout = GaussianOutput("output.log")
final_mol = gout.final_structure
energy = gout.final_energy
frequencies = gout.frequencies
```

### 灯

经典分子动力学。

```python
from pymatgen.io.lammps.data import LammpsData
from pymatgen.io.lammps.inputs import LammpsInputFile

# Structure to LAMMPS data file
lammps_data = LammpsData.from_structure(struct)
lammps_data.write_file("data.lammps")

# LAMMPS input script
lammps_input = LammpsInputFile.from_file("in.lammps")
```

### 量子浓缩咖啡

```python
from pymatgen.io.pwscf import PWInput, PWOutput

# Input
pwin = PWInput(
    struct,
    control={"calculation": "scf"},
    system={"ecutwfc": 50, "ecutrho": 400},
    electrons={"conv_thr": 1e-8}
)
pwin.write_file("pw.in")

# Output
pwout = PWOutput("pw.out")
final_structure = pwout.final_structure
energy = pwout.final_energy
```

### 阿比尼特

```python
from pymatgen.io.abinit import AbinitInput

abin = AbinitInput(struct, pseudos)
abin.set_vars(ecut=10, nband=10)
abin.write("abinit.in")
```

### CP2K

```python
from pymatgen.io.cp2k.inputs import Cp2kInput
from pymatgen.io.cp2k.outputs import Cp2kOutput

# Input
cp2k_input = Cp2kInput.from_file("cp2k.inp")

# Output
cp2k_output = Cp2kOutput("cp2k.out")
```

### FEFF（XAS/XANES）

```python
from pymatgen.io.feff import FeffInput

feff_input = FeffInput(struct, absorbing_atom="Fe")
feff_input.write_file("feff.inp")
```

### LMTO（斯图加特 TB-LMTO-ASA）

```python
from pymatgen.io.lmto import LMTOCtrl

ctrl = LMTOCtrl.from_file("CTRL")
ctrl.structure
```

### Q-化学

```python
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput

# Input
qc_input = QCInput(
    mol,
    rem={"method": "B3LYP", "basis": "6-31G*", "job_type": "opt"}
)
qc_input.write_file("mol.qin")

# Output
qc_output = QCOutput("mol.qout")
```

### 令人兴奋

```python
from pymatgen.io.exciting import ExcitingInput

exc_input = ExcitingInput(struct)
exc_input.write_file("input.xml")
```

### ATAT（合金理论自动化工具包）

```python
from pymatgen.io.atat import Mcsqs

mcsqs = Mcsqs(struct)
mcsqs.write_input(".")
```

## 特殊用途格式

### 盗版

```python
from pymatgen.io.phonopy import get_phonopy_structure, get_pmg_structure

# Convert to phonopy structure
phonopy_struct = get_phonopy_structure(struct)

# Convert from phonopy
struct = get_pmg_structure(phonopy_struct)
```

### ASE（原子模拟环境）

```python
from pymatgen.io.ase import AseAtomsAdaptor

adaptor = AseAtomsAdaptor()

# Pymatgen to ASE
atoms = adaptor.get_atoms(struct)

# ASE to Pymatgen
struct = adaptor.get_structure(atoms)
```

### Zeo++（多孔材料）

```python
from pymatgen.io.zeopp import get_voronoi_nodes, get_high_accuracy_voronoi_nodes

# Analyze pore structure
vor_nodes = get_voronoi_nodes(struct)
```

### BabelMolAdaptor (OpenBabel)

```python
from pymatgen.io.babel import BabelMolAdaptor

adaptor = BabelMolAdaptor(mol)

# Convert to different formats
pdb_str = adaptor.pdbstring
sdf_str = adaptor.write_file("mol.sdf", file_format="sdf")

# Generate 3D coordinates
adaptor.add_hydrogen()
adaptor.make3d()
```

## 炼金术和转换 I/O

### 变形结构

跟踪其转变历史的结构。

```python
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.transformations.standard_transformations import (
    SupercellTransformation,
    SubstitutionTransformation
)

# Create transformed structure
ts = TransformedStructure(struct, [])
ts.append_transformation(SupercellTransformation([[2,0,0],[0,2,0],[0,0,2]]))
ts.append_transformation(SubstitutionTransformation({"Fe": "Mn"}))

# Write with history
ts.write_vasp_input("./calc_dir")

# Read from SNL (Structure Notebook Language)
ts = TransformedStructure.from_snl(snl)
```

## 批量操作

### CifTransmuter

处理多个 CIF 文件。

```python
from pymatgen.alchemy.transmuters import CifTransmuter

transmuter = CifTransmuter.from_filenames(
    ["structure1.cif", "structure2.cif"],
    [SupercellTransformation([[2,0,0],[0,2,0],[0,0,2]])]
)

# Write all structures
transmuter.write_vasp_input("./batch_calc")
```

### PoscarTransmuter

POSCAR 文件类似。

```python
from pymatgen.alchemy.transmuters import PoscarTransmuter

transmuter = PoscarTransmuter.from_filenames(
    ["POSCAR1", "POSCAR2"],
    [transformation1, transformation2]
)
```

## 最佳实践

1. **自动格式检测**：尽可能使用 `from_file()` 和 `to()` 方法
2. **错误处理**：始终将文件 I/O 包装在 try- except 块中
3. **特定于格式的解析器**：使用专门的解析器（例如，`Vasprun`）进行详细的输出分析
4. **输入集**：预配置输入集优于手动参数指定
5. **序列化**：使用JSON/YAML进行长期存储和版本控制
6. **批处理**：使用 Transmuters 将转换应用于多个结构

## 支持的格式摘要

### 结构格式：
CIF、POSCAR/CONTCAR、XYZ、PDB、XSF、PWMAT、Res、CSSR、JSON、YAML

### 电子结构代码：
VASP、Gaussian、LAMMPS、Quantum ESPRESSO、ABINIT、CP2K、FEFF、Q-Chem、LMTO、Exciting、NWChem、AIMS、晶体数据格式

### 分子格式：
XYZ、PDB、MOL、SDF、PQR，通过 OpenBabel（许多其他格式）

### 特殊用途：
Phonopy、ASE、Zeo++、Lobster、BoltzTraP