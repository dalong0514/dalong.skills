<!-- 此文件由机器翻译自 advanced_features.md -->

# 高级功能

## 自定义强制

### 强制类型

FluidSim 支持多种强制机制来维持湍流或驱动特定的动态。

#### 时间相关随机强迫

最常见的持续湍流：

```python
params.forcing.enable = True
params.forcing.type = "tcrandom"
params.forcing.nkmin_forcing = 2  # minimum forced wavenumber
params.forcing.nkmax_forcing = 5  # maximum forced wavenumber
params.forcing.forcing_rate = 1.0  # energy injection rate
params.forcing.tcrandom_time_correlation = 1.0  # correlation time
```

#### 比例强制

保持特定的能量分布：

<<<代码块_1>>>

#### 脚本中的自定义强制

直接在启动脚本中定义强制：

<<<代码块_2>>>

## 自定义初始条件

### 脚本内初始化

完全控制初始字段：

<<<代码块_3>>>

### 层初始化（分层流）

设置密度层：

<<<代码块_4>>>

## 使用 MPI 进行并行计算

### 运行 MPI 模拟

安装 MPI 支持：
<<<代码块_5>>>

使用 MPI 运行：
<<<代码块_6>>>

FluidSim 自动检测 MPI 并分配计算。

### MPI 特定参数

```python
# No special parameters needed
# FluidSim handles domain decomposition automatically

# For very large 3D simulations
params.oper.nx = 512
params.oper.ny = 512
params.oper.nz = 512

# Run with: mpirun -np 64 python script.py
```

### 使用 MPI 输出

输出文件是从 0 级处理器写入的。分析脚本对于串行运行和 MPI 运行的工作方式相同。

## 参数研究

### 运行多个模拟

生成并运行多个参数组合的脚本：

```python
from fluidsim.solvers.ns2d.solver import Simul
import numpy as np

# Parameter ranges
viscosities = [1e-3, 5e-4, 1e-4, 5e-5]
resolutions = [128, 256, 512]

for nu in viscosities:
    for nx in resolutions:
        params = Simul.create_default_params()

        # Configure simulation
        params.oper.nx = params.oper.ny = nx
        params.nu_2 = nu
        params.time_stepping.t_end = 10.0

        # Unique output directory
        params.output.sub_directory = f"nu{nu}_nx{nx}"

        # Run simulation
        sim = Simul(params)
        sim.time_stepping.start()
```

### 集群提交

向集群提交多个作业：

```python
from fluiddyn.clusters.legi import Calcul8 as Cluster

cluster = Cluster()

for nu in viscosities:
    for nx in resolutions:
        script_content = f"""
from fluidsim.solvers.ns2d.solver import Simul

params = Simul.create_default_params()
params.oper.nx = params.oper.ny = {nx}
params.nu_2 = {nu}
params.time_stepping.t_end = 10.0
params.output.sub_directory = "nu{nu}_nx{nx}"

sim = Simul(params)
sim.time_stepping.start()
"""

        with open(f"job_nu{nu}_nx{nx}.py", "w") as f:
            f.write(script_content)

        cluster.submit_script(
            f"job_nu{nu}_nx{nx}.py",
            name_run=f"sim_nu{nu}_nx{nx}",
            nb_nodes=1,
            nb_cores_per_node=24,
            walltime="12:00:00"
        )
```

### 分析参数研究

```python
import os
import pandas as pd
from fluidsim import load_sim_for_plot
import matplotlib.pyplot as plt

results = []

# Collect data from all simulations
for sim_dir in os.listdir("simulations"):
    sim_path = f"simulations/{sim_dir}"
    if not os.path.isdir(sim_path):
        continue

    try:
        sim = load_sim_for_plot(sim_path)

        # Extract parameters
        nu = sim.params.nu_2
        nx = sim.params.oper.nx

        # Extract results
        df = sim.output.spatial_means.load()
        final_energy = df["E"].iloc[-1]
        mean_energy = df["E"].mean()

        results.append({
            "nu": nu,
            "nx": nx,
            "final_energy": final_energy,
            "mean_energy": mean_energy
        })
    except Exception as e:
        print(f"Error loading {sim_dir}: {e}")

# Analyze results
results_df = pd.DataFrame(results)

# Plot results
plt.figure(figsize=(10, 6))
for nx in results_df["nx"].unique():
    subset = results_df[results_df["nx"] == nx]
    plt.plot(subset["nu"], subset["mean_energy"],
             marker="o", label=f"nx={nx}")

plt.xlabel("Viscosity")
plt.ylabel("Mean Energy")
plt.xscale("log")
plt.legend()
plt.savefig("parametric_study_results.png")
```

## 自定义求解器

### 扩展现有求解器

通过继承现有的求解器来创建一个新的求解器：

```python
from fluidsim.solvers.ns2d.solver import Simul as SimulNS2D
from fluidsim.base.setofvariables import SetOfVariables

class SimulCustom(SimulNS2D):
    """Custom solver with additional physics"""

    @staticmethod
    def _complete_params_with_default(params):
        """Add custom parameters"""
        SimulNS2D._complete_params_with_default(params)
        params._set_child("custom", {"param1": 0.0})

    def __init__(self, params):
        super().__init__(params)
        # Custom initialization

    def tendencies_nonlin(self, state_spect=None):
        """Override to add custom tendencies"""
        tendencies = super().tendencies_nonlin(state_spect)

        # Add custom terms
        # tendencies.vx_fft += custom_term_vx
        # tendencies.vy_fft += custom_term_vy

        return tendencies
```

使用自定义求解器：
```python
params = SimulCustom.create_default_params()
# Configure params...
sim = SimulCustom(params)
sim.time_stepping.start()
```

## 在线可视化

仿真期间显示字段：

```python
params.output.ONLINE_PLOT_OK = True
params.output.periods_plot.phys_fields = 1.0  # plot every 1.0 time units
params.output.phys_fields.field_to_plot = "vorticity"

sim = Simul(params)
sim.time_stepping.start()
```

绘图在执行过程中实时出现。

## 检查点并重新启动

### 自动检查点

```python
params.output.periods_save.phys_fields = 1.0  # save every 1.0 time units
```

模拟过程中会自动保存字段。

### 手动检查点

```python
# During simulation
sim.output.phys_fields.save()
```

### 从检查点重新启动

```python
params = Simul.create_default_params()
params.init_fields.type = "from_file"
params.init_fields.from_file.path = "simulation_dir/state_phys_t5.000.h5"
params.time_stepping.t_end = 20.0  # extend simulation

sim = Simul(params)
sim.time_stepping.start()
```

## 内存和性能优化

### 减少内存使用

```python
# Disable unnecessary outputs
params.output.periods_save.spectra = 0  # disable spectra saving
params.output.periods_save.spect_energy_budg = 0  # disable energy budget

# Reduce spatial field saves
params.output.periods_save.phys_fields = 10.0  # save less frequently
```

### 优化 FFT 性能

```python
import os

# Select FFT library
os.environ["FLUIDSIM_TYPE_FFT2D"] = "fft2d.with_fftw"
os.environ["FLUIDSIM_TYPE_FFT3D"] = "fft3d.with_fftw"

# Or use MKL if available
# os.environ["FLUIDSIM_TYPE_FFT2D"] = "fft2d.with_mkl"
```

### 时间步长优化

```python
# Use adaptive time stepping
params.time_stepping.USE_CFL = True
params.time_stepping.CFL = 0.8  # slightly larger CFL for faster runs

# Use efficient time scheme
params.time_stepping.type_time_scheme = "RK4"  # 4th order Runge-Kutta
```