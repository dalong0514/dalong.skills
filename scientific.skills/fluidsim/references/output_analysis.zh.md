<!-- 此文件由机器翻译自 output_analysis.md -->

# 输出与分析

## 输出类型

FluidSim 在模拟过程中自动保存多种类型的输出。

### 物理场

**文件格式**：HDF5 (`.h5`)

**位置**：`simulation_dir/state_phys_t*.h5`

**内容**：特定时间的速度、涡度和其他物理空间场

**访问**：
```python
sim.output.phys_fields.plot()
sim.output.phys_fields.plot("vorticity")
sim.output.phys_fields.plot("vx")
sim.output.phys_fields.plot("div")  # check divergence

# Save manually
sim.output.phys_fields.save()

# Get data
vorticity = sim.state.state_phys.get_var("rot")
```

### 空间手段

**文件格式**：文本文件 (`.txt`)

**位置**：`simulation_dir/spatial_means.txt`

**内容**：体积平均量与时间的关系（能量、熵等）

**访问**：
<<<代码块_1>>>

### 光谱

**文件格式**：HDF5 (`.h5`)

**位置**：`simulation_dir/spectra_*.h5`

**内容**：能量和熵谱与波数

**访问**：
<<<代码块_2>>>

### 光谱能量预算

**文件格式**：HDF5 (`.h5`)

**位置**：`simulation_dir/spect_energy_budg_*.h5`

**内容**：尺度之间的能量传递

**访问**：
<<<代码块_3>>>

## 后处理

### 加载模拟进行分析

#### 快速加载（只读）

<<<代码块_4>>>

使用它进行快速可视化和分析。不初始化完整模拟状态。

#### 完整状态加载

<<<代码块_5>>>

### 可视化工具

#### 内置绘图

FluidSim 通过 matplotlib 提供基本绘图：

<<<代码块_6>>>

#### 高级可视化

对于出版质量或 3D 可视化：

**ParaView**：直接打开`.h5`文件
```bash
paraview simulation_dir/state_phys_t*.h5
```

**VisIt**：类似于大型数据集的 ParaView

**自定义Python**：
```python
import h5py
import matplotlib.pyplot as plt

# Load field manually
with h5py.File("state_phys_t10.000.h5", "r") as f:
    vx = f["state_phys"]["vx"][:]
    vy = f["state_phys"]["vy"][:]

# Custom plotting
plt.contourf(vx)
plt.show()
```

## 分析实例

### 能源进化

```python
from fluidsim import load_sim_for_plot
import matplotlib.pyplot as plt

sim = load_sim_for_plot("simulation_dir")
df = sim.output.spatial_means.load()

plt.figure()
plt.plot(df["t"], df["E"], label="Kinetic Energy")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.show()
```

### 光谱分析

```python
sim = load_sim_for_plot("simulation_dir")

# Plot energy spectrum
sim.output.spectra.plot1d(tmin=5.0, tmax=10.0)  # average over time range

# Get spectral data
k, E_k = sim.output.spectra.load1d_mean(tmin=5.0, tmax=10.0)

# Check for power law
import numpy as np
log_k = np.log(k)
log_E = np.log(E_k)
# fit power law in inertial range
```

### 参数研究分析

使用不同参数运行多个模拟时：

```python
import os
import pandas as pd
from fluidsim import load_sim_for_plot

# Collect results from multiple simulations
results = []
for sim_dir in os.listdir("simulations"):
    if not os.path.isdir(f"simulations/{sim_dir}"):
        continue

    sim = load_sim_for_plot(f"simulations/{sim_dir}")

    # Extract key metrics
    df = sim.output.spatial_means.load()
    final_energy = df["E"].iloc[-1]

    # Get parameters
    nu = sim.params.nu_2

    results.append({
        "nu": nu,
        "final_energy": final_energy,
        "sim_dir": sim_dir
    })

# Analyze results
results_df = pd.DataFrame(results)
results_df.plot(x="nu", y="final_energy", logx=True)
```

### 字段操作

```python
sim = load_sim_for_plot("simulation_dir")

# Load specific time
sim.output.phys_fields.set_of_phys_files.update_times()
times = sim.output.phys_fields.set_of_phys_files.times

# Load field at specific time
field_file = sim.output.phys_fields.get_field_to_plot(time=5.0)
vorticity = field_file.get_var("rot")

# Compute derived quantities
import numpy as np
vorticity_rms = np.sqrt(np.mean(vorticity**2))
vorticity_max = np.max(np.abs(vorticity))
```

## 输出目录结构

```
simulation_dir/
├── params_simul.xml         # Simulation parameters
├── stdout.txt               # Standard output log
├── state_phys_t*.h5         # Physical fields at different times
├── spatial_means.txt        # Time series of spatial averages
├── spectra_*.h5            # Spectral data
├── spect_energy_budg_*.h5  # Energy budget data
└── info_solver.txt         # Solver information
```

## 性能监控

```python
# During simulation, check progress
sim.output.print_stdout.complete_timestep()

# After simulation, review performance
sim.output.print_stdout.plot_deltat()  # plot time step evolution
sim.output.print_stdout.plot_clock_times()  # plot computation time
```

## 数据导出

将 Fluidsim 输出转换为其他格式：

```python
import h5py
import numpy as np

# Export to numpy array
with h5py.File("state_phys_t10.000.h5", "r") as f:
    vx = f["state_phys"]["vx"][:]
    np.save("vx.npy", vx)

# Export to CSV
df = sim.output.spatial_means.load()
df.to_csv("spatial_means.csv", index=False)
```