<!-- 此文件由机器翻译自 material-handling.md -->

# PyLabRobot 中的物料搬运设备

## 概述

PyLabRobot 与材料处理设备集成，包括加热摇床、培养箱、离心机和泵。这些设备可实现基本液体处理之外的环境控制、样品制备和自动化工作流程。

## 加热摇床

### 汉密尔顿 HeaterShaker

Hamilton HeaterShaker 为微孔板提供温度控制和轨道振荡。

#### 设置

```python
from pylabrobot.heating_shaking import HeaterShaker
from pylabrobot.heating_shaking.hamilton import HamiltonHeaterShakerBackend

# Create heater shaker
hs = HeaterShaker(
    name="heater_shaker_1",
    backend=HamiltonHeaterShakerBackend(),
    size_x=156.0,
    size_y=  156.0,
    size_z=18.0
)

await hs.setup()
```

#### 操作

**温度控制：**

<<<代码块_1>>>

**震动控制：**

<<<代码块_2>>>

**板块操作：**

<<<代码块_3>>>

#### 与液体处理器集成

<<<代码块_4>>>

#### 多个加热摇床

HamiltonHeaterShakerBackend 处理多个单元：

<<<代码块_5>>>

### Inheco ThermoShake

Inheco ThermoShake 提供温度控制和摇动。

#### 设置

<<<代码块_6>>>

#### 操作

类似于 Hamilton HeaterShaker：

```python
# Temperature control
await hs.set_temperature(37)
temp = await hs.get_temperature()

# Shaking control
await hs.set_shake_rate(300)

# Plate locking
await hs.lock_plate()
await hs.unlock_plate()
```

## 孵化器

### Inheco孵化器

PyLabRobot 支持各种 Inheco 培养箱型号，用于温控板存储。

#### 支持的型号

- Inheco 单板培养箱
- Inheco 多板培养箱
- 其他 Inheco 温度控制器

#### 设置

```python
from pylabrobot.temperature_control import TemperatureController
from pylabrobot.temperature_control.inheco import InhecoBackend

# Create incubator
incubator = TemperatureController(
    name="incubator",
    backend=InhecoBackend(),
    size_x=156.0,
    size_y=156.0,
    size_z=50.0
)

await incubator.setup()
```

#### 操作

```python
# Set temperature
await incubator.set_temperature(37)

# Get temperature
temp = await incubator.get_temperature()
print(f"Incubator temperature: {temp}°C")

# Turn off
await incubator.set_temperature(None)
```

### Thermo Fisher Cytomat 培养箱

Cytomat 培养箱提供带有温度和 CO2 控制的自动化板存储。

#### 设置

```python
from pylabrobot.incubation import Incubator
from pylabrobot.incubation.cytomat_backend import CytomatBackend

incubator = Incubator(
    name="cytomat",
    backend=CytomatBackend()
)

await incubator.setup()
```

#### 操作

```python
# Store plate
await incubator.store_plate(plate_id="plate_001", position=1)

# Retrieve plate
await incubator.retrieve_plate(position=1)

# Set environmental conditions
await incubator.set_temperature(37)
await incubator.set_co2(5.0)  # 5% CO2
```

## 离心机

### 安捷伦 VSpin

Agilent VSpin 是一款用于板处理的真空辅助离心机。

#### 设置

```python
from pylabrobot.centrifuge import Centrifuge
from pylabrobot.centrifuge.vspin import VSpinBackend

centrifuge = Centrifuge(
    name="vspin",
    backend=VSpinBackend()
)

await centrifuge.setup()
```

#### 操作

**门控制：**

```python
# Open door
await centrifuge.open_door()

# Close door
await centrifuge.close_door()

# Lock door
await centrifuge.lock_door()

# Unlock door
await centrifuge.unlock_door()
```

**铲斗定位：**

```python
# Move bucket to loading position
await centrifuge.move_bucket_to_loading()

# Move bucket to home position
await centrifuge.move_bucket_to_home()
```

**旋转：**

```python
# Run centrifuge
await centrifuge.spin(
    speed=2000,      # RPM
    duration=300     # seconds
)

# Stop spinning
await centrifuge.stop_spin()
```

#### 集成示例

```python
async def centrifuge_workflow():
    """Complete centrifugation workflow"""

    lh = LiquidHandler(backend=STAR(), deck=STARLetDeck())
    centrifuge = Centrifuge(name="vspin", backend=VSpinBackend())

    await lh.setup()
    await centrifuge.setup()

    try:
        # Prepare samples
        await lh.pick_up_tips(tip_rack["A1:H1"])
        await lh.transfer(samples["A1:H12"], plate["A1:H12"], vols=200)
        await lh.drop_tips()

        # Load into centrifuge
        print("Move plate to centrifuge")
        await centrifuge.open_door()
        await centrifuge.move_bucket_to_loading()
        input("Press Enter when plate is loaded...")

        await centrifuge.move_bucket_to_home()
        await centrifuge.close_door()
        await centrifuge.lock_door()

        # Centrifuge
        await centrifuge.spin(speed=2000, duration=300)

        # Unload
        await centrifuge.unlock_door()
        await centrifuge.open_door()
        await centrifuge.move_bucket_to_loading()
        input("Press Enter when plate is removed...")

        await centrifuge.move_bucket_to_home()
        await centrifuge.close_door()

    finally:
        await lh.stop()
        await centrifuge.stop()
```

## 泵

### 科尔·帕默 Masterflex

PyLabRobot 支持 Cole Parmer Masterflex 蠕动泵进行流体传输。

#### 设置

```python
from pylabrobot.pumps import Pump
from pylabrobot.pumps.cole_parmer import ColeParmerMasterflexBackend

pump = Pump(
    name="masterflex",
    backend=ColeParmerMasterflexBackend()
)

await pump.setup()
```

#### 操作

**运行泵：**

```python
# Run for duration
await pump.run_for_duration(
    duration=10,      # seconds
    speed=50          # % of maximum
)

# Run continuously
await pump.start(speed=50)

# Stop pump
await pump.stop()
```

**基于体积的泵送：**

```python
# Pump specific volume (requires calibration)
await pump.pump_volume(
    volume=10,        # mL
    speed=50          # % of maximum
)
```

#### 校准

```python
# Calibrate pump for volume accuracy
# (requires known volume measurement)
await pump.run_for_duration(duration=60, speed=50)
actual_volume = 25.3  # mL measured

pump.calibrate(duration=60, speed=50, volume=actual_volume)
```

### Agrowtek 泵阵列

支持 Agrowtek 泵阵列，可同时进行多个流体传输。

#### 设置

```python
from pylabrobot.pumps import PumpArray
from pylabrobot.pumps.agrowtek import AgrowtekBackend

pump_array = PumpArray(
    name="agrowtek",
    backend=AgrowtekBackend(),
    num_pumps=8
)

await pump_array.setup()
```

#### 操作

```python
# Run specific pump
await pump_array.run_pump(
    pump_number=1,
    duration=10,
    speed=50
)

# Run multiple pumps simultaneously
await pump_array.run_pumps(
    pump_numbers=[1, 2, 3],
    duration=10,
    speed=50
)
```

## 多设备协议

### 复杂工作流程示例

```python
async def complex_workflow():
    """Multi-device automated workflow"""

    # Initialize all devices
    lh = LiquidHandler(backend=STAR(), deck=STARLetDeck())
    hs = HeaterShaker(name="hs", backend=HamiltonHeaterShakerBackend())
    centrifuge = Centrifuge(name="vspin", backend=VSpinBackend())
    pump = Pump(name="pump", backend=ColeParmerMasterflexBackend())

    await lh.setup()
    await hs.setup()
    await centrifuge.setup()
    await pump.setup()

    try:
        # 1. Sample preparation
        await lh.pick_up_tips(tip_rack["A1:H1"])
        await lh.transfer(samples["A1:H12"], plate["A1:H12"], vols=100)
        await lh.drop_tips()

        # 2. Add reagent via pump
        await pump.pump_volume(volume=50, speed=50)

        # 3. Mix on heater shaker
        await hs.lock_plate()
        await hs.set_temperature(37)
        await hs.set_shake_rate(300)
        await asyncio.sleep(600)  # 10 min incubation
        await hs.set_shake_rate(0)
        await hs.set_temperature(None)
        await hs.unlock_plate()

        # 4. Centrifuge
        await centrifuge.open_door()
        # (load plate)
        await centrifuge.close_door()
        await centrifuge.spin(speed=2000, duration=180)
        await centrifuge.open_door()
        # (unload plate)

        # 5. Transfer supernatant
        await lh.pick_up_tips(tip_rack["A2:H2"])
        await lh.transfer(
            plate["A1:H12"],
            output_plate["A1:H12"],
            vols=80
        )
        await lh.drop_tips()

    finally:
        await lh.stop()
        await hs.stop()
        await centrifuge.stop()
        await pump.stop()
```

## 最佳实践

1. **设备初始化**：在协议启动时设置所有设备
2. **顺序操作**：物料搬运通常需要顺序步骤
3. **安全**：在手动处理板材之前务必解锁/打开门
4. **温度平衡**：让设备有时间达到温度
5. **错误处理**：使用try/finally优雅地处理设备错误
6. **状态验证**：操作前检查设备状态
7. **计时**：考虑设备特定的延迟（加热、离心）
8. **维护**：遵循制造商维护计划
9. **校准**：定期校准泵和温度控制器
10. **文档**：记录所有设备设置和参数

## 常见模式

### 温控孵化

```python
async def incubate_with_shaking(
    plate,
    temperature: float,
    shake_rate: int,
    duration: int
):
    """Incubate plate with temperature and shaking"""

    hs = HeaterShaker(name="hs", backend=HamiltonHeaterShakerBackend())
    await hs.setup()

    try:
        # Assign plate to heater shaker
        hs.assign_child_resource(plate, location=(0, 0, 0))

        # Start incubation
        await hs.lock_plate()
        await hs.set_temperature(temperature)
        await hs.set_shake_rate(shake_rate)

        # Wait
        await asyncio.sleep(duration)

        # Stop
        await hs.set_shake_rate(0)
        await hs.set_temperature(None)
        await hs.unlock_plate()

    finally:
        await hs.stop()

# Use in protocol
await incubate_with_shaking(
    plate=assay_plate,
    temperature=37,
    shake_rate=300,
    duration=600  # 10 minutes
)
```

### 自动化板材加工

```python
async def process_plates(plate_list: list):
    """Process multiple plates through workflow"""

    lh = LiquidHandler(backend=STAR(), deck=STARLetDeck())
    hs = HeaterShaker(name="hs", backend=HamiltonHeaterShakerBackend())

    await lh.setup()
    await hs.setup()

    try:
        for i, plate in enumerate(plate_list):
            print(f"Processing plate {i+1}/{len(plate_list)}")

            # Transfer samples
            await lh.pick_up_tips(tip_rack[f"A{i+1}:H{i+1}"])
            await lh.transfer(
                source[f"A{i+1}:H{i+1}"],
                plate["A1:H1"],
                vols=100
            )
            await lh.drop_tips()

            # Incubate
            hs.assign_child_resource(plate, location=(0, 0, 0))
            await hs.lock_plate()
            await hs.set_temperature(37)
            await hs.set_shake_rate(300)
            await asyncio.sleep(300)  # 5 min
            await hs.set_shake_rate(0)
            await hs.set_temperature(None)
            await hs.unlock_plate()
            hs.unassign_child_resource(plate)

    finally:
        await lh.stop()
        await hs.stop()
```

## 其他资源

- 物料搬运文档：https://docs.pylabrobot.org/user_guide/01_material-handling/
- 加热器摇床：https://docs.pylabrobot.org/user_guide/01_material-handling/heating_shaking/
- API 参考：https://docs.pylabrobot.org/api/
- 支持的设备：https://docs.pylabrobot.org/user_guide/machines.html