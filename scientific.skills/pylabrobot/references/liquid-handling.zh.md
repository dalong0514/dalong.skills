<!-- 此文件由机器翻译自 liquid-handling.md -->

# 使用 PyLabRobot 处理液体

## 概述

液体处理模块 (`pylabrobot.liquid_handling`) 提供用于控制液体处理机器人的统一接口。 `LiquidHandler` 类充当所有移液操作的主接口，通过后端抽象跨不同的硬件平台工作。

## 基本设置

### 初始化液体处理器

```python
from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends import STAR
from pylabrobot.resources import STARLetDeck

# Create liquid handler with STAR backend
lh = LiquidHandler(backend=STAR(), deck=STARLetDeck())
await lh.setup()

# When done
await lh.stop()
```

### 在后端之间切换

通过交换后端来更改机器人，无需重写协议：

<<<代码块_1>>>

## 核心运营

### 小费管理

拾取和滴下吸头对于液体处理操作至关重要：

<<<代码块_2>>>

**小费跟踪**：启用自动小费跟踪以监控小费使用情况：

<<<代码块_3>>>

### 吸取液体

从孔或容器中抽取液体：

<<<代码块_4>>>

### 分配液体

将液体分配到孔或容器中：

<<<代码块_5>>>

### 转移液体

转移将吸液和分配结合在一次操作中：

<<<代码块_6>>>

## 先进技术

### 系列稀释

跨板行或列创建连续稀释：

```python
# 2-fold serial dilution
source_vols = [100, 50, 50, 50, 50, 50, 50, 50]
dest_vols = [0, 50, 50, 50, 50, 50, 50, 50]

# Add diluent first
await lh.pick_up_tips(tip_rack["A1"])
await lh.transfer(
    source=buffer["A1"],
    dest=plate["A2:A8"],
    vols=50
)
await lh.drop_tips()

# Perform serial dilution
await lh.pick_up_tips(tip_rack["A2"])
for i in range(7):
    await lh.aspirate(plate[f"A{i+1}"], vols=50)
    await lh.dispense(plate[f"A{i+2}"], vols=50)
    # Mix
    await lh.aspirate(plate[f"A{i+2}"], vols=50)
    await lh.dispense(plate[f"A{i+2}"], vols=50)
await lh.drop_tips()
```

### 板复制

将整个板布局复制到另一个板：

```python
# Setup tips
await lh.pick_up_tips(tip_rack["A1:H1"])

# Replicate 96-well plate (12 columns)
for col in range(1, 13):
    await lh.transfer(
        source=source_plate[f"A{col}:H{col}"],
        dest=dest_plate[f"A{col}:H{col}"],
        vols=100
    )

await lh.drop_tips()
```

### 多通道移液

同时使用多个通道进行并行操作：

```python
# 8-channel transfer (entire row)
await lh.pick_up_tips(tip_rack["A1:H1"])
await lh.transfer(
    source=source_plate["A1:H1"],
    dest=dest_plate["A1:H1"],
    vols=100
)
await lh.drop_tips()

# Process entire plate with 8-channel
for col in range(1, 13):
    await lh.pick_up_tips(tip_rack[f"A{col}:H{col}"])
    await lh.transfer(
        source=source_plate[f"A{col}:H{col}"],
        dest=dest_plate[f"A{col}:H{col}"],
        vols=100
    )
    await lh.drop_tips()
```

### 混合液体

通过反复吸取和分配来混合液体：

```python
# Mix by aspiration/dispensing
await lh.pick_up_tips(tip_rack["A1"])

# Mix 5 times
for _ in range(5):
    await lh.aspirate(plate["A1"], vols=80)
    await lh.dispense(plate["A1"], vols=80)

await lh.drop_tips()
```

## 音量追踪

自动跟踪孔中的液体体积：

```python
from pylabrobot.resources import set_volume_tracking

# Enable volume tracking globally
set_volume_tracking(True)

# Set initial volumes
plate["A1"].tracker.set_liquids([(None, 200)])  # 200 µL

# After aspirating 100 µL
await lh.aspirate(plate["A1"], vols=100)
print(plate["A1"].tracker.get_volume())  # 100 µL

# Check remaining volume
remaining = plate["A1"].tracker.get_volume()
```

## 液体课程

定义最佳移液的液体特性：

```python
# Liquid classes control aspiration/dispense parameters
from pylabrobot.liquid_handling import LiquidClass

# Create custom liquid class
water = LiquidClass(
    name="Water",
    aspiration_flow_rate=100,
    dispense_flow_rate=150,
    aspiration_mix_flow_rate=100,
    dispense_mix_flow_rate=100,
    air_transport_retract_dist=10
)

# Use with operations
await lh.aspirate(
    plate["A1"],
    vols=100,
    liquid_class=water
)
```

## 错误处理

处理液体处理操作中的错误：

```python
try:
    await lh.setup()
    await lh.pick_up_tips(tip_rack["A1"])
    await lh.transfer(source["A1"], dest["A1"], vols=100)
    await lh.drop_tips()
except Exception as e:
    print(f"Error during liquid handling: {e}")
    # Attempt to drop tips if holding them
    try:
        await lh.drop_tips()
    except:
        pass
finally:
    await lh.stop()
```

## 最佳实践

1. **始终设置和停止**：在操作之前调用 `await lh.setup()` 并在完成后调用 `await lh.stop()`
2. **启用跟踪**：使用尖端跟踪和音量跟踪进行准确的状态管理
3. **吸头管理**：吸液前务必拿起吸头，完成后将其放下
4. **流速**：根据液体粘度和容器类型调整流速
5. **液体高度**：设置适当的吸取/分配高度以避免飞溅
6. **错误处理**：使用 try/finally 块来确保正确的清理
7. **模拟测试**：在硬件上运行之前使用 ChatterboxBackend 测试协议
8. **体积限制**：遵守吸头体积限制和孔容量
9. **混合**：在分配粘性液体后或在精确度至关重要时进行混合
10. **文档**：记录液体类别和自定义参数以实现可重复性

## 常见模式

### 完整的液体处理协议

```python
from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends import STAR
from pylabrobot.resources import STARLetDeck, TIP_CAR_480_A00, Cos_96_DW_1mL
from pylabrobot.resources import set_tip_tracking, set_volume_tracking

# Enable tracking
set_tip_tracking(True)
set_volume_tracking(True)

# Initialize
lh = LiquidHandler(backend=STAR(), deck=STARLetDeck())
await lh.setup()

try:
    # Define resources
    tip_rack = TIP_CAR_480_A00(name="tips")
    source = Cos_96_DW_1mL(name="source")
    dest = Cos_96_DW_1mL(name="dest")

    # Assign to deck
    lh.deck.assign_child_resource(tip_rack, rails=1)
    lh.deck.assign_child_resource(source, rails=10)
    lh.deck.assign_child_resource(dest, rails=15)

    # Set initial volumes
    for well in source.children:
        well.tracker.set_liquids([(None, 200)])

    # Execute protocol
    await lh.pick_up_tips(tip_rack["A1:H1"])
    await lh.transfer(
        source=source["A1:H12"],
        dest=dest["A1:H12"],
        vols=100
    )
    await lh.drop_tips()

finally:
    await lh.stop()
```

## 硬件特定注释

### 汉密尔顿之星

- 支持完整的液体处理能力
- 使用USB连接进行通信
- 直接执行固件命令
- 支持 CO-RE（压缩 O 形环扩张）尖端

### Opentrons OT-2

- 需要 IP 地址才能进行网络连接
- 使用HTTP API进行通信
- 仅限于 8 通道和单通道移液器
- 与 STAR 相比，甲板布局更简单

### 帝肯 EVO

- 正在进行的工作支持
- 与 Hamilton STAR 类似的功能
- 检查文档中当前的兼容性状态

## 其他资源

- 官方液体处理指南：https://docs.pylabrobot.org/user_guide/basic.html
- API 参考：https://docs.pylabrobot.org/api/pylabrobot.liquid_handling.html
- 协议示例：https://github.com/PyLabRobot/pylabrobot/tree/main/examples