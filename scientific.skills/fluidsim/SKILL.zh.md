<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：流体模拟
描述：使用 Python 进行计算流体动力学模拟的框架。在运行流体动力学模拟（包括纳维-斯托克斯方程 (2D/3D)、浅水方程、分层流）或分析湍流、涡动力学或地球物理流时使用。提供具有 FFT、HPC 支持和全面输出分析的伪谱方法。
---

# 流体模拟

## 概述

FluidSim 是一个面向对象的 Python 框架，用于高性能计算流体动力学 (CFD) 模拟。它使用伪谱方法和 FFT 提供周期域方程的求解器，提供与 Fortran/C++ 相当的性能，同时保持 Python 的易用性。

**主要优势**：
- 多个求解器：2D/3D Navier-Stokes、浅水、分层流
- 高性能：Pythran/Transonic 编译、MPI 并行化
- 完整的工作流程：参数配置、仿真执行、输出分析
- 交互式分析：基于Python的后处理和可视化

## 核心能力

### 1. 安装和设置

使用 uv 和适当的功能标志安装 Fluidsim：

```bash
# Basic installation
uv uv pip install fluidsim

# With FFT support (required for most solvers)
uv uv pip install "fluidsim[fft]"

# With MPI for parallel computing
uv uv pip install "fluidsim[fft,mpi]"
```

设置输出目录的环境变量（可选）：

<<<代码块_1>>>

无需 API 密钥或身份验证。

有关完整的安装说明和环境配置，请参阅`references/installation.md`。

### 2. 运行模拟

标准工作流程包括五个步骤：

**第 1 步**：导入求解器
<<<代码块_2>>>

**第2步**：创建并配置参数
<<<代码块_3>>>

**步骤 3**：实例化模拟
<<<代码块_4>>>

**第4步**：执行
<<<代码块_5>>>

**第 5 步**：分析结果
<<<代码块_6>>>

有关完整示例、重新启动模拟和集群部署，请参阅 `references/simulation_workflow.md`。

### 3.可用的求解器

根据物理问题选择求解器：

**2D Navier-Stokes** (`ns2d`)：2D 湍流、涡流动力学
```python
from fluidsim.solvers.ns2d.solver import Simul
```

**3D Navier-Stokes** (`ns3d`)：3D 湍流，真实流动
```python
from fluidsim.solvers.ns3d.solver import Simul
```

**分层流** (`ns2d.strat`, `ns3d.strat`)：海洋/大气流
```python
from fluidsim.solvers.ns2d.strat.solver import Simul
params.N = 1.0  # Brunt-Väisälä frequency
```

**浅水** (`sw1l`)：地球物理流、旋转系统
```python
from fluidsim.solvers.sw1l.solver import Simul
params.f = 1.0  # Coriolis parameter
```

请参阅 `references/solvers.md` 以获取完整的求解器列表和选择指南。

### 4.参数配置

参数按层次结构组织并通过点表示法访问：

**域名和分辨率**：
```python
params.oper.nx = 256  # grid points
params.oper.Lx = 2 * pi  # domain size
```

**物理参数**：
```python
params.nu_2 = 1e-3  # viscosity
params.nu_4 = 0     # hyperviscosity (optional)
```

**时间步进**：
```python
params.time_stepping.t_end = 10.0
params.time_stepping.USE_CFL = True  # adaptive time step
params.time_stepping.CFL = 0.5
```

**初始条件**：
```python
params.init_fields.type = "noise"  # or "dipole", "vortex", "from_file", "in_script"
```

**输出设置**：
```python
params.output.periods_save.phys_fields = 1.0  # save every 1.0 time units
params.output.periods_save.spectra = 0.5
params.output.periods_save.spatial_means = 0.1
```

对于拼写错误，Parameters 对象会引发 `AttributeError`，以防止静默配置错误。

请参阅 `references/parameters.md` 以获取全面的参数文档。

### 5. 输出与分析

FluidSim 产生多种输出类型，在模拟过程中自动保存：

**物理场**：HDF5 格式的速度、涡度
```python
sim.output.phys_fields.plot("vorticity")
sim.output.phys_fields.plot("vx")
```

**空间手段**：体积平均数量的时间序列
```python
sim.output.spatial_means.plot()
```

**光谱**：能量和熵光谱
```python
sim.output.spectra.plot1d()
sim.output.spectra.plot2d()
```

**加载之前的模拟**：
```python
from fluidsim import load_sim_for_plot
sim = load_sim_for_plot("simulation_dir")
sim.output.phys_fields.plot()
```

**高级可视化**：在 ParaView 或 VisIt 中打开 `.h5` 文件进行 3D 可视化。

请参阅 `references/output_analysis.md` 了解详细的分析工作流程、参数研究分析和数据导出。

### 6. 高级功能

**自定义强迫**：维持湍流或驱动特定动态
```python
params.forcing.enable = True
params.forcing.type = "tcrandom"  # time-correlated random forcing
params.forcing.forcing_rate = 1.0
```

**自定义初始条件**：在脚本中定义字段
```python
params.init_fields.type = "in_script"
sim = Simul(params)
X, Y = sim.oper.get_XY_loc()
vx = sim.state.state_phys.get_var("vx")
vx[:] = sin(X) * cos(Y)
sim.time_stepping.start()
```

**MPI 并行化**：在多个处理器上运行
```bash
mpirun -np 8 python simulation_script.py
```

**参数研究**：使用不同参数运行多个模拟
```python
for nu in [1e-3, 5e-4, 1e-4]:
    params = Simul.create_default_params()
    params.nu_2 = nu
    params.output.sub_directory = f"nu{nu}"
    sim = Simul(params)
    sim.time_stepping.start()
```

请参阅 `references/advanced_features.md` 了解强制类型、自定义求解器、集群提交和性能优化。

## 常见用例

### 二维湍流研究

```python
from fluidsim.solvers.ns2d.solver import Simul
from math import pi

params = Simul.create_default_params()
params.oper.nx = params.oper.ny = 512
params.oper.Lx = params.oper.Ly = 2 * pi
params.nu_2 = 1e-4
params.time_stepping.t_end = 50.0
params.time_stepping.USE_CFL = True
params.init_fields.type = "noise"
params.output.periods_save.phys_fields = 5.0
params.output.periods_save.spectra = 1.0

sim = Simul(params)
sim.time_stepping.start()

# Analyze energy cascade
sim.output.spectra.plot1d(tmin=30.0, tmax=50.0)
```

### 分层流模拟

```python
from fluidsim.solvers.ns2d.strat.solver import Simul

params = Simul.create_default_params()
params.oper.nx = params.oper.ny = 256
params.N = 2.0  # stratification strength
params.nu_2 = 5e-4
params.time_stepping.t_end = 20.0

# Initialize with dense layer
params.init_fields.type = "in_script"
sim = Simul(params)
X, Y = sim.oper.get_XY_loc()
b = sim.state.state_phys.get_var("b")
b[:] = exp(-((X - 3.14)**2 + (Y - 3.14)**2) / 0.5)
sim.state.statephys_from_statespect()

sim.time_stepping.start()
sim.output.phys_fields.plot("b")
```

### 使用 MPI 进行高分辨率 3D 模拟

```python
from fluidsim.solvers.ns3d.solver import Simul

params = Simul.create_default_params()
params.oper.nx = params.oper.ny = params.oper.nz = 512
params.nu_2 = 1e-5
params.time_stepping.t_end = 10.0
params.init_fields.type = "noise"

sim = Simul(params)
sim.time_stepping.start()
```

运行：
```bash
mpirun -np 64 python script.py
```

### 泰勒-格林涡验证

```python
from fluidsim.solvers.ns2d.solver import Simul
import numpy as np
from math import pi

params = Simul.create_default_params()
params.oper.nx = params.oper.ny = 128
params.oper.Lx = params.oper.Ly = 2 * pi
params.nu_2 = 1e-3
params.time_stepping.t_end = 10.0
params.init_fields.type = "in_script"

sim = Simul(params)
X, Y = sim.oper.get_XY_loc()
vx = sim.state.state_phys.get_var("vx")
vy = sim.state.state_phys.get_var("vy")
vx[:] = np.sin(X) * np.cos(Y)
vy[:] = -np.cos(X) * np.sin(Y)
sim.state.statephys_from_statespect()

sim.time_stepping.start()

# Validate energy decay
df = sim.output.spatial_means.load()
# Compare with analytical solution
```

## 快速参考

**导入求解器**：`from fluidsim.solvers.ns2d.solver import Simul`

**创建参数**：`params = Simul.create_default_params()`

**设置分辨率**：`params.oper.nx = params.oper.ny = 256`

**设置粘度**：`params.nu_2 = 1e-3`

**设置结束时间**：`params.time_stepping.t_end = 10.0`

**运行模拟**：`sim = Simul(params); sim.time_stepping.start()`

**绘图结果**：`sim.output.phys_fields.plot("vorticity")`
**负载模拟**：`sim = load_sim_for_plot("path/to/sim")`

## 资源

**文档**：https://fluidsim.readthedocs.io/

**参考文件**：
- `references/installation.md`：完整的安装说明
- `references/solvers.md`：可用求解器和选择指南
- `references/simulation_workflow.md`：详细的工作流程示例
- `references/parameters.md`：综合参数文档
- `references/output_analysis.md`：输出类型和分析方法
- `references/advanced_features.md`：强制、MPI、参数研究、自定义求解器