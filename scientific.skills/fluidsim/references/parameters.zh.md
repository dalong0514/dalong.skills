<!-- 此文件由机器翻译自 parameters.md -->

# 参数配置

## 参数对象

`Parameters` 对象是分层的并组织成逻辑组。使用点表示法访问：

```python
params = Simul.create_default_params()
params.group.subgroup.parameter = value
```

## 关键参数组

### 运算符 (`params.oper`)

定义域和分辨率：

<<<代码块_1>>>

**分辨率指南**：使用 2 的幂以获得最佳 FFT 性能（128、256、512、1024 等）

### 物理参数

#### 粘度

<<<代码块_2>>>

高阶粘度（`nu_4`、`nu_8`）可抑制高波数而不影响大尺度。

#### 分层（分层求解器）

<<<代码块_3>>>

#### 旋转（浅水）

<<<代码块_4>>>

### 时间步进 (`params.time_stepping`)

<<<代码块_5>>>

**推荐**：将 `USE_CFL=True` 与 `CFL=0.5` 一起使用以实现自适应时间步进。

### 初始字段 (`params.init_fields`)

<<<代码块_6>>>

**可用类型**：
- `"noise"`：随机噪声
- `"dipole"`：涡旋偶极子
- `"vortex"`：单涡流
- `"taylor_green"`：泰勒-格林涡旋
- `"from_file"`：从文件加载
- `"in_script"`：在脚本中定义

#### 来自文件

```python
params.init_fields.type = "from_file"
params.init_fields.from_file.path = "path/to/state_file.h5"
```

#### 在脚本中

```python
params.init_fields.type = "in_script"

# Define initialization after creating sim
sim = Simul(params)

# Access state fields
vx = sim.state.state_phys.get_var("vx")
vy = sim.state.state_phys.get_var("vy")

# Set fields
X, Y = sim.oper.get_XY_loc()
vx[:] = np.sin(X) * np.cos(Y)
vy[:] = -np.cos(X) * np.sin(Y)

# Run simulation
sim.time_stepping.start()
```

### 输出设置 (`params.output`)

#### 输出目录

```python
params.output.sub_directory = "my_simulation"
```

在 `$FLUIDSIM_PATH` 或当前目录中创建的目录。

#### 保存句号

```python
params.output.periods_save.phys_fields = 1.0  # save fields every 1.0 time units
params.output.periods_save.spectra = 0.5      # save spectra
params.output.periods_save.spatial_means = 0.1  # save spatial averages
params.output.periods_save.spect_energy_budg = 0.5  # spectral energy budget
```

设置为 `0` 以禁用特定输出类型。

#### 打印控制

```python
params.output.periods_print.print_stdout = 0.5  # print status every 0.5 time units
```

#### 在线绘图

```python
params.output.periods_plot.phys_fields = 2.0  # plot every 2.0 time units

# Must also enable the output module
params.output.ONLINE_PLOT_OK = True
params.output.phys_fields.field_to_plot = "vorticity"  # or "vx", "vy", etc.
```

### 强制 (`params.forcing`)

添加强迫项以维持能量：

```python
params.forcing.enable = True
params.forcing.type = "tcrandom"  # time-correlated random forcing

# Forcing parameters
params.forcing.nkmax_forcing = 5  # maximum forced wavenumber
params.forcing.nkmin_forcing = 2  # minimum forced wavenumber
params.forcing.forcing_rate = 1.0  # energy injection rate
```

**常见的强迫类型**：
- `"tcrandom"`：时间相关的随机强迫
- `"proportional"`：比例强制（保持特定频谱）
- `"in_script"`：脚本中定义的自定义强制

## 参数安全

当访问不存在的参数时，Parameters 对象会引发 `AttributeError`：

```python
params.nu_2 = 1e-3  # OK
params.nu2 = 1e-3   # ERROR: AttributeError
```

这可以防止在基于文本的配置文件中被默默忽略的拼写错误。

## 查看所有参数

```python
# Print all parameters
params._print_as_xml()

# Get as dictionary
param_dict = params._make_dict()
```

## 保存参数配置

参数与模拟输出一起自动保存：

```python
params._save_as_xml("simulation_params.xml")
params._save_as_json("simulation_params.json")
```