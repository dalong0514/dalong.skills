<!-- 此文件由机器翻译自 chemistry_molecular_formats.md -->

# 化学和分子文件格式参考

本参考涵盖计算化学、化学信息学、分子建模和相关领域中常用的文件格式。

## 结构文件格式

### .pdb - 蛋白质数据库
**描述：**生物大分子3D结构的标准格式
**典型数据：** 原子坐标、残基信息、二级结构、晶体结构数据
**用例：** 蛋白质结构分析、分子可视化、对接研究
**Python 库：**
- `Biopython`：`Bio.PDB`
- `MDAnalysis`：`MDAnalysis.Universe('file.pdb')`
- `PyMOL`：`pymol.cmd.load('file.pdb')`
- `ProDy`: `prody.parsePDB('file.pdb')`
**EDA 方法：**
- 结构验证（键长、角度、碰撞）
- 二级结构分析
- B 因子分布
- 缺失残留物/原子检测
- 用于验证的 Ramachandran 图
- 表面积和体积计算

### .cif - 晶体信息文件
**描述：** 晶体学信息的结构化数据格式
**典型数据：**晶胞参数、原子坐标、对称操作、实验数据
**用例：** 晶体结构测定、结构生物学、材料科学
**Python 库：**
- `gemmi`: `gemmi.cif.read_file('file.cif')`
- `PyCifRW`：`CifFile.ReadCif('file.cif')`
- `Biopython`：`Bio.PDB.MMCIFParser()`
**EDA 方法：**
- 数据完整性检查
- 分辨率和质量指标
- 晶胞参数分析
- 对称群验证
- 原子位移参数
- R 因素和验证指标

### .mol - MDL Molfile
**描述：** MDL/Accelrys 的化学结构文件格式
**典型数据：** 2D/3D 坐标、原子类型、键序、电荷
**用例：** 化学数据库存储、化学信息学、药物设计
**Python 库：**
- `RDKit`：`Chem.MolFromMolFile('file.mol')`
- `Open Babel`：`pybel.readfile('mol', 'file.mol')`
- `ChemoPy`：用于描述符计算
**EDA 方法：**
- 分子性质计算（MW、logP、TPSA）
- 官能团分析
- 环系统检测
- 立体化学验证
- 2D/3D坐标一致性
- 化合价和电荷验证

### .mol2 - Tripos Mol2
**描述：** 具有原子类型的完整 3D 分子结构格式
**典型数据：** 坐标、SYBYL 原子类型、键类型、电荷、子结构
**用例：** 分子对接、QSAR 研究、药物发现
**Python 库：**
- `RDKit`：`Chem.MolFromMol2File('file.mol2')`
- `Open Babel`：`pybel.readfile('mol2', 'file.mol2')`
- `MDAnalysis`：可以解析 mol2 拓扑
**EDA 方法：**
- 原子类型分布
- 部分电荷分析
- 债券类型统计
- 子结构识别
- 构象分析
- 能量最小化状态检查

### .sdf - 结构数据文件
**描述：** 具有关联数据的多结构文件格式
**典型数据：** 具有属性/注释的多个分子结构
**用例：** 化学数据库、虚拟筛选、化合物库
**Python 库：**
- `RDKit`：`Chem.SDMolSupplier('file.sdf')`
- `Open Babel`：`pybel.readfile('sdf', 'file.sdf')`
- `PandasTools` (RDKit)：用于 DataFrame 集成
**EDA 方法：**
- 数据集大小和多样性指标
- 性能分布分析（MW、logP等）
- 结构多样性（谷本相似性）
- 缺失数据评估
- 属性中的异常值检测
- 支架分析

### .xyz - XYZ 坐标
**说明：** 简单笛卡尔坐标格式
**典型数据：** 原子类型和 3D 坐标
**用例：** 量子化学、几何优化、分子动力学
**Python 库：**
- `ASE`：`ase.io.read('file.xyz')`
- `Open Babel`：`pybel.readfile('xyz', 'file.xyz')`
- `cclib`：用于使用 xyz 解析 QM 输出
**EDA 方法：**
- 几何分析（键长、角度、二面角）
- 质心计算
- 转动惯量
- 分子大小指标
- 坐标验证
- 对称性检测

### .smi / .smiles - SMILES 字符串
**描述：**化学结构的行符号
**典型数据：** 分子结构的文本表示
**用例：** 化学数据库、文献挖掘、数据交换
**Python 库：**
- `RDKit`：`Chem.MolFromSmiles(smiles)`
- `Open Babel`：可以解析 SMILES
- `DeepChem`：用于 SMILES 上的 ML
**EDA 方法：**
- SMILES语法验证
- SMILES 的描述符计算
- 指纹生成
- 子结构搜索
- 互变异构体计数
- 立体异构体处理

### .pdbqt - AutoDock PDBQT
**说明：** 修改PDB格式，用于AutoDock对接
**典型数据：** 坐标、部分电荷、对接原子类型
**用例：** 分子对接、虚拟筛选
**Python 库：**
- `Meeko`：用于 PDBQT 准备
- `Open Babel`：可以读取 PDBQT
- `ProDy`：有限的 PDBQT 支持
**EDA 方法：**
- 电荷分布分析
- 可旋转键识别
- 原子类型验证
- 协调质量检查
- 氢放置验证
- 扭转定义分析

### .mae - Maestro 格式
**描述：** 薛定谔专有的分子结构格式
**典型数据：** 薛定谔套件的结构、属性、注释
**用例：** 药物发现、使用薛定谔工具进行分子建模
**Python 库：**
- `schrodinger.structure`：需要安装薛定谔
- 用于基本阅读的自定义解析器
**EDA 方法：**
- 属性提取和分析
- 结构质量指标
- 符合者分析
- 对接得分分布
- 配体效率指标

### .gro - GROMACS 坐标文件
**描述：** GROMACS MD 模拟的分子结构文件
**典型数据：** 原子位置、速度、框向量
**用例：** 分子动力学模拟、GROMACS 工作流程
**Python 库：**
- `MDAnalysis`：`Universe('file.gro')`
- `MDTraj`：`mdtraj.load_gro('file.gro')`
- `GromacsWrapper`：用于 GROMACS 集成
**EDA 方法：**
- 系统构成分析
- 盒子尺寸验证
- 原子位置分布
- 速度分布（如果存在）
- 密度计算
- 溶剂化分析

## 计算化学输出格式

### .log - 高斯日志文件
**描述：** 高斯量子化学计算的输出
**典型数据：** 能量、几何形状、频率、轨道、人口
**用例：** QM 计算、几何优化、频率分析
**Python 库：**
- `cclib`：`cclib.io.ccread('file.log')`
- `GaussianRunPack`：对于高斯工作流程
- 带有正则表达式的自定义解析器
**EDA 方法：**
- 收敛性分析
- 能量剖面提取
- 振动频率分析
- 轨道能级
- 人口分析（Mulliken，NBO）
- 热化学数据提取

### .out - 量子化学输出
**描述：** 来自各种 QM 包的通用输出文件
**典型数据：** 计算结果、能量、性质
**用例：**跨不同软件的 QM 计算
**Python 库：**
- `cclib`：QM 输出的通用解析器
- `ASE`：可以读取一些输出格式
**EDA 方法：**
- 软件特定的解析
- 收敛标准检查
- 能量和梯度趋势
- 基础集和方法验证
- 计算成本分析

### .wfn / .wfx - 波函数文件
**描述：** 用于量子化学分析的波函数数据
**典型数据：** 分子轨道、基组、密度矩阵
**用例：** 电子密度分析、QTAIM 分析
**Python 库：**
- `Multiwfn`：通过 Python 的接口
- `Horton`：用于波函数分析
- 针对特定格式的自定义解析器
**EDA 方法：**
- 轨道人口分析
- 电子密度分布
- 临界点分析（QTAIM）
- 分子轨道可视化
- 粘合分析

### .fchk - 高斯格式检查点
**描述：** Gaussian 格式化的检查点文件
**典型数据：** 完整的波函数数据、结果、几何形状
**用例：** 后处理高斯计算
**Python 库：**
- `cclib`：可以解析fchk文件
- `GaussView` Python API（如果可用）
- 自定义解析器
**EDA 方法：**
- 波函数质量评估
- 财产提取
- 基组信息
- 梯度和Hessian分析
- 自然轨道分析

### .cube - 高斯立方体文件
**描述：** 3D 网格上的体积数据
**典型数据：** 电子密度、分子轨道、网格上的 ESP
**用例：**体积属性的可视化
**Python 库：**
- `cclib`：`cclib.io.ccread('file.cube')`
- `ase.io`：`ase.io.read('file.cube')`
- `pyquante`：用于多维数据集文件操作
**EDA 方法：**
- 网格尺寸和间距分析
- 价值分布统计
- 等值面值测定
- 体积积分
- 不同立方体之间的比较
## 分子动力学格式

### .dcd - 二元轨迹
**描述：** 二元轨迹格式（CHARMM、NAMD）
**典型数据：**原子坐标的时间序列
**用例：** MD 轨迹分析
**Python 库：**
- `MDAnalysis`：`Universe(topology, 'traj.dcd')`
- `MDTraj`：`mdtraj.load_dcd('traj.dcd', top='topology.pdb')`
- `PyTraj`（琥珀色）：有限支持
**EDA 方法：**
- RMSD/RMSF 分析
- 轨迹长度和帧数
- 坐标范围和漂移
- 定期边界处理
- 文件完整性检查
- 时间步验证

### .xtc - 压缩轨迹
**说明：** GROMACS 压缩轨迹格式
**典型数据：** MD 模拟的压缩坐标
**用例：** 节省空间的 MD 轨迹存储
**Python 库：**
- `MDAnalysis`：`Universe(topology, 'traj.xtc')`
- `MDTraj`：`mdtraj.load_xtc('traj.xtc', top='topology.pdb')`
**EDA 方法：**
- 压缩比评估
- 精度损失评估
- 随着时间的推移 RMSD
- 结构稳定性指标
- 采样频率分析

### .trr - GROMACS 轨迹
**描述：** 全精度 GROMACS 轨迹
**典型数据：** 纵向坐标、速度、力
**用例：** 高精度 MD 分析
**Python 库：**
- `MDAnalysis`：完全支持
- `MDTraj`：可以读取 trr 文件
- `GromacsWrapper`
**EDA 方法：**
- 完整的系统动力学分析
- 能量守恒检查（速度）
- 力分析
- 温度和压力验证
- 系统平衡评估

### .nc / .netcdf - 琥珀色 NetCDF 轨迹
**描述：** 网络通用数据形式轨迹
**典型数据：** MD 坐标、速度、力
**用例：** Amber MD 模拟、大轨迹存储
**Python 库：**
- `MDAnalysis`：NetCDF 支持
- `PyTraj`：原生琥珀分析
- `netCDF4`：低级访问
**EDA 方法：**
- 元数据提取
- 轨迹统计
- 时间序列分析
- 副本交换分析
- 多维度数据提取

### .top - GROMACS 拓扑
**描述：** GROMACS 的分子拓扑
**典型数据：** 原子类型、键、角度、力场参数
**用例：** MD 模拟设置和分析
**Python 库：**
- `ParmEd`：`parmed.load_file('system.top')`
- `MDAnalysis`：可以解析拓扑
- 针对特定字段的自定义解析器
**EDA 方法：**
- 力场参数验证
- 系统组成
- 键/角/二面角分布
- 电荷中性检查
- 分子类型计数

### .psf - 蛋白质结构文件 (CHARMM)
**描述：** CHARMM/NAMD 的拓扑文件
**典型数据：** Atom 连接性、类型、费用
**用例：** CHARMM/NAMD MD 模拟
**Python 库：**
- `MDAnalysis`：本机 PSF 支持
- `ParmEd`：可以读取 PSF 文件
**EDA 方法：**
- 拓扑验证
- 连接性分析
- 电荷分布
- 原子类型统计
- 细分分析

### .prmtop - 琥珀色参数/拓扑
**描述：** Amber拓扑和参数文件
**典型数据：**系统拓扑、力场参数
**用例：** Amber MD 模拟
**Python 库：**
- `ParmEd`：`parmed.load_file('system.prmtop')`
- `PyTraj`：本机 Amber 支持
**EDA 方法：**
- 力场完整性
- 参数验证
- 系统规模和组成
- 定期盒子信息
- 创建原子掩模以进行分析

### .inpcrd / .rst7 - 琥珀色坐标
**描述：** 琥珀色坐标/重启文件
**典型数据：** 原子坐标、速度、盒子信息
**用例：** Amber MD 的起始坐标
**Python 库：**
- `ParmEd`：与 prmtop 一起使用
- `PyTraj`：琥珀色坐标读取
**EDA 方法：**
- 坐标有效性
- 系统初始化检查
- 框向量验证
- 速度分布（如果重新启动）
- 能源最小化状态

## 光谱和分析数据

### .jcamp / .jdx - JCAMP-DX
**描述：**原子和分子物理数据交换联合委员会
**典型数据：** 光谱数据（IR、NMR、MS、UV-Vis）
**用例：** 光谱数据交换和归档
**Python 库：**
- `jcamp`：`jcamp.jcamp_reader('file.jdx')`
- `nmrglue`：对于 NMR JCAMP 文件
- 针对特定子类型的自定义解析器
**EDA 方法：**
- 峰值检测和分析
- 基线校正评估
- 信噪比计算
- 光谱范围验证
- 整合分析
- 与参考光谱的比较
### .mzML - 质谱标记语言
**描述：** 质谱数据的标准 XML 格式
**典型数据：** MS/MS 谱图、色谱图、元数据
**用例：** 蛋白质组学、代谢组学、质谱工作流程
**Python 库：**
- `pymzml`：`pymzml.run.Reader('file.mzML')`
- `pyteomics`：`pyteomics.mzml.read('file.mzML')`
- `MSFileReader` 包装器
**EDA 方法：**
- 扫描数量和类型
- MS等级分布
- 保留时间范围
- m/z 范围和分辨率
- 峰值强度分布
- 数据完整性
- 质量控制指标

### .mzXML - 质谱 XML
**描述：** MS 数据的开放 XML 格式
**典型数据：** 质谱、保留时间、峰列表
**用例：** 传统 MS 数据、代谢组学
**Python 库：**
- `pymzml`：可以读取 mzXML
- `pyteomics.mzxml`
- `lxml` 用于直接 XML 解析
**EDA 方法：**
- 类似于mzML
- 版本兼容性检查
- 转换质量评估
- 峰值选取验证

### .raw - 供应商原始数据
**描述：** 专有仪器数据文件（Thermo、Bruker 等）
**典型数据：** 原始仪器信号，未处理的数据
**用例：** 直接仪器数据访问
**Python 库：**
- `pymsfilereader`：对于 Thermo RAW 文件
- `ThermoRawFileParser`：CLI 包装器
- 供应商特定的 API（Thermo、Bruker Compass）
**EDA 方法：**
- 仪器法提取
- 原始信号质量
- 校准状态
- 扫描功能分析
- 色谱质量指标

### .d - 安捷伦数据目录
**说明：** 安捷伦的数据文件夹结构
**典型数据：** LC-MS、GC-MS 数据和元数据
**用例：** 安捷伦仪器数据处理
**Python 库：**
- `agilent-reader`：社区工具
- `Chemstation` Python 集成
- 自定义目录解析
**EDA 方法：**
- 目录结构验证
- 方法参数提取
- 信号文件完整性
- 校准曲线分析
- 序列信息提取

### .fid - 无核磁共振感应衰变
**描述：**原始 NMR 时域数据
**典型数据：** 时域 NMR 信号
**用例：** NMR 处理和分析
**Python 库：**
- `nmrglue`：`nmrglue.bruker.read_fid('fid')`
- `nmrstarlib`：对于 NMR-STAR 文件
**EDA 方法：**
- 信号衰减分析
- 噪音水平评估
- 采集参数验证
- 变迹函数选择
- 零填充优化
- 定相参数估计

### .ft - NMR 频域数据
**描述：** 处理后的核磁共振谱
**典型数据：**频域NMR数据
**用例：** NMR 分析和解释
**Python 库：**
- `nmrglue`：全面的 NMR 支持
- `pyNMR`：用于处理
**EDA 方法：**
- 峰值选取和整合
- 化学位移校准
- 多重性分析
- 耦合常数提取
- 频谱质量指标
- 参考化合物鉴定

### .spc - 光谱文件
**描述：**热银河光谱格式
**典型数据：** 红外、拉曼、紫外-可见光谱
**用例：** 来自各种仪器的光谱数据
**Python 库：**
- `spc`：`spc.File('file.spc')`
- 二进制格式的自定义解析器
**EDA 方法：**
- 光谱分辨率
- 波长/波数范围
- 基线表征
- 峰值识别
- 导数光谱计算

## 化学数据库格式

### .inchi - 国际化学标识符
**描述：**化学物质的文本标识符
**典型数据：** 层状化学结构表示
**用例：** 化学数据库密钥、结构搜索
**Python 库：**
- `RDKit`：`Chem.MolFromInchi(inchi)`
- `Open Babel`：InChI 转换
**EDA 方法：**
- InChI验证
- 层分析
- 立体化学验证
- InChI 密钥生成
- 结构往返验证

### .cdx / .cdxml - ChemDraw 交换
**说明：** ChemDraw 绘图文件格式
**典型数据：** 带注释的二维化学结构
**用例：** 化学绘图、出版物数据
**Python 库：**
- `RDKit`：可以导入一些CDXML
- `Open Babel`：有限支持
- `ChemDraw` Python API（商业）
**EDA 方法：**
- 结构提取
- 注释保存
- 风格一致性
- 2D坐标验证

### .cml - 化学标记语言
**描述：** 基于XML的化学结构格式
**典型数据：** 化学结构、反应、性质
**用例：** 语义化学数据表示
**Python 库：**
- `RDKit`：CML 支持
- `Open Babel`：良好的 CML 支持
- `lxml`：用于 XML 解析
**EDA 方法：**
- XML 模式验证
- 命名空间处理
- 财产提取
- 反应方案分析
- 元数据完整性

### .rxn - MDL 反应文件
**说明：**化学反应结构文件
**典型数据：**反应物、产物、反应箭头
**用例：** 反应数据库、合成规划
**Python 库：**
- `RDKit`：`Chem.ReactionFromRxnFile('file.rxn')`
- `Open Babel`：反应支持
**EDA 方法：**
- 反应平衡验证
- 原子图谱分析
- 试剂鉴定
- 立体化学变化
- 反应分类

### .rdf - 反应数据文件
**说明：**多反应文件格式
**典型数据：** 多重反应数据
**用例：** 反应数据库
**Python 库：**
- `RDKit`：RDF 读取功能
- 自定义解析器
**EDA 方法：**
- 反应收率统计
- 状态分析
- 成功率模式
- 试剂频率分析

## 计算输出和数据

### .hdf5 / .h5 - 分层数据格式
**描述：**科学数据数组的容器
**典型数据：** 大型数组、元数据、分层组织
**用例：** 大型数据集存储、计算结果
**Python 库：**
- `h5py`：`h5py.File('file.h5', 'r')`
- `pytables`：高级 HDF5 接口
- `pandas`：可以读取 HDF5
**EDA 方法：**
- 数据集结构探索
- 数组形状和数据类型分析
- 元数据提取
- 内存高效的数据采样
- 块优化分析
- 压缩比评估

### .pkl / .pickle - Python Pickle
**描述：** 序列化的Python对象
**典型数据：** 任何 Python 对象（分子、数据帧、模型）
**用例：** 中间数据存储、模型持久化
**Python 库：**
- `pickle`：内置序列化
- `joblib`：增强大型阵列的酸洗
- `dill`：扩展pickle支持
**EDA 方法：**
- 对象类型检查
- 规模和复杂性分析
- 版本兼容性检查
- 安全验证（可信来源）
- 反序列化测试

### .npy / .npz - NumPy 数组
**说明：** NumPy 数组二进制格式
**典型数据：** 数值数组（坐标、特征、矩阵）
**用例：** 快速数值数据 I/O
**Python 库：**
- `numpy`：`np.load('file.npy')`
- 大文件的直接内存映射
**EDA 方法：**
- 阵列形状和尺寸
- 数据类型和精度
- 统计摘要（平均值、标准差、范围）
- 缺失值检测
- 异常值识别
- 内存占用分析

### .mat - MATLAB 数据文件
**描述：** MATLAB 工作区数据
**典型数据：** MATLAB 中的数组、结构
**用例：** MATLAB-Python 数据交换
**Python 库：**
- `scipy.io`：`scipy.io.loadmat('file.mat')`
- `h5py`：适用于 v7.3 MAT 文件
**EDA 方法：**
- 变量提取和类型
- 数组维数分析
- 构造领域勘探
- MATLAB版本兼容性
- 数据类型转换验证

### .csv - 逗号分隔值
**描述：** 文本格式的表格数据
**典型数据：** 化学性质、实验数据、描述符
**用例：** 数据交换、分析、机器学习
**Python 库：**
- `pandas`：`pd.read_csv('file.csv')`
- `csv`：内置模块
- `polars`：快速读取 CSV
**EDA 方法：**
- 数据类型推断
- 缺失值模式
- 统计摘要
- 相关性分析
- 分布可视化
- 异常值检测

### .json - JavaScript 对象表示法
**描述：** 结构化文本数据格式
**典型数据：** 化学性质、元数据、API 响应
**用例：** 数据交换、配置、Web API
**Python 库：**
- `json`：内置 JSON 支持
- `pandas`：`pd.read_json()`
- `ujson`：更快的 JSON 解析
**EDA 方法：**
- 模式验证
- 嵌套深度分析
- 键值分布
- 数据类型一致性
- 数组长度统计

### .parquet - Apache Parquet
**说明：** 列式存储格式
**典型数据：** 高效的大型表格数据集
**用例：** 大数据，高效的柱状分析
**Python 库：**
- `pandas`：`pd.read_parquet('file.parquet')`
- `pyarrow`：直接镶木地板访问
- `fastparquet`：替代实现
**EDA 方法：**
- 来自元数据的列统计信息
- 分区分析
- 压缩效率
- 行组结构
- 大文件快速采样
- 模式演变跟踪