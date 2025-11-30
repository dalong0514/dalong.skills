<!-- 此文件由机器翻译自 analytical-equipment.md -->

# PyLabRobot 中的分析设备

## 概述

PyLabRobot 与分析设备集成，包括读板器、秤和其他测量设备。这允许将液体处理与分析测量相结合的自动化工作流程。

## 读板器

### BMG CLARIOstar（增强版）

BMG Labtech CLARIOstar 和 CLARIOstar Plus 是用于测量吸光度、发光和荧光的酶标仪。

#### 硬件设置

**物理连接：**
1. IEC C13 电源线连接至主电源
2. 连接计算机的 USB-B 电缆（设备端带有安全螺丝）
3.可选：用于板堆叠单元的RS-232端口

**通讯：**
- 在固件级别通过 FTDI/USB-A 进行串行连接
- 跨平台支持（Windows、macOS、Linux）

#### 软件设置

```python
from pylabrobot.plate_reading import PlateReader
from pylabrobot.plate_reading.clario_star_backend import CLARIOstarBackend

# Create backend
backend = CLARIOstarBackend()

# Initialize plate reader
pr = PlateReader(
    name="CLARIOstar",
    backend=backend,
    size_x=0.0,    # Physical dimensions not critical for plate readers
    size_y=0.0,
    size_z=0.0
)

# Setup (initializes device)
await pr.setup()

# When done
await pr.stop()
```

#### 基本操作

**开幕和闭幕：**

<<<代码块_1>>>

**温度控制：**

<<<代码块_2>>>

**读数测量值：**

<<<代码块_3>>>

#### 数据格式

读板器方法返回数组数据：

<<<代码块_4>>>

#### 与液体处理器集成

将板读数与液体处理相结合：

<<<代码块_5>>>

#### 高级功能

**发展状况：**

一些 CLARIOstar 功能正在开发中：
- 光谱扫描
- 注射针控制
- 详细的测量参数配置
- 特定的阅读模式

检查当前文档以获取最新功能支持。

#### 最佳实践

1. **温度控制**：加热慢，尽早设定温度
2. **板装载**：关闭前确保板正确就位
3. **测量选择**：为您的测定选择合适的波长
4. **数据验证**：检查测量质量和预期范围
5. **错误处理**：处理超时和通信错误
6. **维护**：按照制造商指南保持光学器件清洁

#### 示例：完整的车牌读取工作流程

<<<代码块_6>>>

## 音阶

### 梅特勒托利多秤

PyLabRobot 支持梅特勒托利多秤进行质量测量。

#### 设置

```python
from pylabrobot.scales import Scale
from pylabrobot.scales.mettler_toledo_backend import MettlerToledoBackend

# Create scale
scale = Scale(
    name="analytical_scale",
    backend=MettlerToledoBackend()
)

await scale.setup()
```

#### 操作

```python
# Get weight measurement
weight = await scale.get_weight()  # Returns weight in grams
print(f"Weight: {weight} g")

# Tare (zero) the scale
await scale.tare()

# Get multiple measurements
weights = []
for i in range(5):
    w = await scale.get_weight()
    weights.append(w)
    await asyncio.sleep(1)

average_weight = sum(weights) / len(weights)
print(f"Average weight: {average_weight} g")
```

#### 与液体处理器集成

```python
# Weigh samples during protocol
lh = LiquidHandler(backend=STAR(), deck=STARLetDeck())
scale = Scale(name="scale", backend=MettlerToledoBackend())

await lh.setup()
await scale.setup()

try:
    # Tare scale
    await scale.tare()

    # Dispense liquid
    await lh.pick_up_tips(tip_rack["A1"])
    await lh.aspirate(reagent["A1"], vols=1000)

    # (Move to scale position)

    # Dispense and weigh
    await lh.dispense(container, vols=1000)
    weight = await scale.get_weight()

    print(f"Dispensed weight: {weight} g")

    # Calculate actual volume (assuming density = 1 g/mL for water)
    actual_volume = weight * 1000  # Convert g to µL
    print(f"Actual volume: {actual_volume} µL")

    await lh.drop_tips()

finally:
    await lh.stop()
    await scale.stop()
```

## 其他分析设备

### 流式细胞仪

一些流式细胞仪集成正在开发中。检查当前文档以了解支持状态。

### 分光光度计

可能支持其他分光光度计型号。检查当前设备兼容性的文档。

## 多设备工作流程

### 协调多个设备

```python
async def multi_device_workflow():
    """Coordinate liquid handler, plate reader, and scale"""

    # Initialize all devices
    lh = LiquidHandler(backend=STAR(), deck=STARLetDeck())
    pr = PlateReader(name="CLARIOstar", backend=CLARIOstarBackend())
    scale = Scale(name="scale", backend=MettlerToledoBackend())

    await lh.setup()
    await pr.setup()
    await scale.setup()

    try:
        # 1. Weigh reagent
        await scale.tare()
        # (place container on scale)
        reagent_weight = await scale.get_weight()

        # 2. Prepare samples with liquid handler
        await lh.pick_up_tips(tip_rack["A1:H1"])
        await lh.transfer(source["A1:H12"], dest["A1:H12"], vols=100)
        await lh.drop_tips()

        # 3. Read plate
        await pr.open()
        # (load plate)
        await pr.close()
        data = await pr.read_absorbance(wavelength=450)

        return {
            "reagent_weight": reagent_weight,
            "absorbance_data": data
        }

    finally:
        await lh.stop()
        await pr.stop()
        await scale.stop()
```

## 最佳实践

1. **设备初始化**：在协议开始时设置所有设备
2. **错误处理**：优雅地处理通信错误
3. **清理**：始终在所有设备上调用 `stop()`
4. **计时**：考虑设备特定的计时（温度平衡、测量时间）
5. **校准**：遵循制造商校准程序
6. **数据验证**：验证测量结果是否在预期范围内
7. **文档**：记录设备设置和参数
8. **集成测试**：彻底测试多设备工作流程
9. **并发操作**：尽可能使用异步来重叠操作
10. **数据存储**：使用元数据（时间戳、设置）保存原始数据

## 常见模式

### 动态板读数

```python
async def kinetic_reading(num_reads: int, interval: int):
    """Perform kinetic plate reading"""

    pr = PlateReader(name="CLARIOstar", backend=CLARIOstarBackend())
    await pr.setup()

    try:
        await pr.set_temperature(37)
        await pr.open()
        # (load plate)
        await pr.close()

        results = []
        for i in range(num_reads):
            data = await pr.read_absorbance(wavelength=450)
            timestamp = time.time()
            results.append({
                "read_number": i + 1,
                "timestamp": timestamp,
                "data": data
            })

            if i < num_reads - 1:
                await asyncio.sleep(interval)

        return results

    finally:
        await pr.stop()

# Read every 30 seconds for 10 minutes
results = await kinetic_reading(num_reads=20, interval=30)
```

## 其他资源

- 车牌读取文档：https://docs.pylabrobot.org/user_guide/02_analytical/
- BMG CLARIOstar 指南：https://docs.pylabrobot.org/user_guide/02_analytical/plate-reading/bmg-clariostar.html
- API 参考：https://docs.pylabrobot.org/api/pylabrobot.plate_reading.html
- 支持的设备：https://docs.pylabrobot.org/user_guide/machines.html