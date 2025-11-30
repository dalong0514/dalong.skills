<!-- 此文件由机器翻译自 visualization.md -->

# PyLabRobot 中的可视化和模拟

## 概述

PyLabRobot 提供可视化和模拟工具，用于在没有物理硬件的情况下开发、测试和验证实验室协议。可视化工具提供甲板状态的实时 3D 可视化，而模拟后端则支持协议测试和验证。

## 展示台

### 展示台是什么？

PyLabRobot Visualizer 是一个基于浏览器的工具，它：
- 显示甲板布局的 3D 可视化
- 显示实时尖端存在和液体体积
- 适用于模拟和物理机器人
- 提供交互式甲板状态检查
- 启用可视化协议验证

### 启动可视化工具

可视化工具作为 Web 服务器运行并显示在您的浏览器中：

```python
from pylabrobot.visualizer import Visualizer

# Create visualizer
vis = Visualizer()

# Start web server (opens browser automatically)
await vis.start()

# Stop visualizer
await vis.stop()
```

**默认设置：**
- 端口：1234 (http://localhost:1234)
- 启动时自动打开浏览器

### 将液体处理器连接到展示台

<<<代码块_1>>>

### 追踪功能

#### 启用跟踪

要使可视化工具显示吸头和液体，请启用跟踪：

<<<代码块_2>>>

#### 设置初始液体

定义初始液体含量以进行可视化：

<<<代码块_3>>>

#### 可视化提示存在

<<<代码块_4>>>

### 完整的可视化示例

<<<代码块_5>>>

## 牌组布局编辑器

### 使用牌组编辑器

PyLabRobot 包含一个图形甲板布局编辑器：

**特点：**
- 视觉甲板设计界面
- 拖放资源放置
- 编辑初始液态
- 设置提示存在
- 将布局保存/加载为 JSON

**用途：**
- 通过可视化界面访问
- 以图形方式而不是代码创建布局
- 导出为 JSON 以在协议中使用

### 加载牌组布局

<<<代码块_6>>>

## 模拟

### 话匣子后端

ChatterboxBackend 模拟液体处理操作：

**特点：**
- 无需硬件
- 验证协议逻辑
- 跟踪提示和音量
- 支持所有液体处理操作
- 与可视化工具一起使用

**设置：**

```python
from pylabrobot.liquid_handling.backends.simulation import ChatterboxBackend

# Create simulation backend
backend = ChatterboxBackend(
    num_channels=8  # Simulate 8-channel pipette
)

# Use with liquid handler
lh = LiquidHandler(backend=backend, deck=STARLetDeck())
```

### 模拟用例

#### 协议开发

```python
async def develop_protocol():
    """Develop protocol using simulation"""

    # Use simulation for development
    lh = LiquidHandler(
        backend=ChatterboxBackend(),
        deck=STARLetDeck()
    )

    # Connect visualizer
    vis = Visualizer()
    await vis.start()
    lh.visualizer = vis

    await lh.setup()

    try:
        # Develop and test protocol
        await lh.pick_up_tips(tip_rack["A1"])
        await lh.transfer(plate["A1"], plate["A2"], vols=100)
        await lh.drop_tips()

        print("Protocol development complete!")

    finally:
        await lh.stop()
        await vis.stop()
```

#### 协议验证

```python
async def validate_protocol():
    """Validate protocol logic without hardware"""

    set_tip_tracking(True)
    set_volume_tracking(True)

    lh = LiquidHandler(
        backend=ChatterboxBackend(),
        deck=STARLetDeck()
    )
    await lh.setup()

    try:
        # Setup resources
        tip_rack = TIP_CAR_480_A00(name="tips")
        plate = Cos_96_DW_1mL(name="plate")

        lh.deck.assign_child_resource(tip_rack, rails=1)
        lh.deck.assign_child_resource(plate, rails=10)

        # Set initial state
        for well in plate.children:
            well.tracker.set_liquids([(None, 200)])

        # Execute protocol
        await lh.pick_up_tips(tip_rack["A1:H1"])

        # Test different volumes
        test_volumes = [50, 100, 150]
        for i, vol in enumerate(test_volumes):
            await lh.transfer(
                plate[f"A{i+1}:H{i+1}"],
                plate[f"A{i+4}:H{i+4}"],
                vols=vol
            )

        await lh.drop_tips()

        # Validate volumes
        for i, vol in enumerate(test_volumes):
            for row in "ABCDEFGH":
                well = plate[f"{row}{i+4}"]
                actual_vol = well.tracker.get_volume()
                assert actual_vol == vol, f"Volume mismatch in {well.name}"

        print("✓ Protocol validation passed!")

    finally:
        await lh.stop()
```

#### 测试边缘情况

```python
async def test_edge_cases():
    """Test protocol edge cases in simulation"""

    lh = LiquidHandler(
        backend=ChatterboxBackend(),
        deck=STARLetDeck()
    )
    await lh.setup()

    try:
        # Test 1: Empty well aspiration
        try:
            await lh.aspirate(empty_plate["A1"], vols=100)
            print("✗ Should have raised error for empty well")
        except Exception as e:
            print(f"✓ Correctly raised error: {e}")

        # Test 2: Overfilling well
        try:
            await lh.dispense(small_well, vols=1000)  # Too much
            print("✗ Should have raised error for overfilling")
        except Exception as e:
            print(f"✓ Correctly raised error: {e}")

        # Test 3: Tip capacity
        try:
            await lh.aspirate(large_volume_well, vols=2000)  # Exceeds tip capacity
            print("✗ Should have raised error for tip capacity")
        except Exception as e:
            print(f"✓ Correctly raised error: {e}")

    finally:
        await lh.stop()
```

### CI/CD 集成

使用模拟进行自动化测试：

```python
# test_protocols.py
import pytest
from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends.simulation import ChatterboxBackend

@pytest.mark.asyncio
async def test_transfer_protocol():
    """Test liquid transfer protocol"""

    lh = LiquidHandler(
        backend=ChatterboxBackend(),
        deck=STARLetDeck()
    )
    await lh.setup()

    try:
        # Setup
        tip_rack = TIP_CAR_480_A00(name="tips")
        plate = Cos_96_DW_1mL(name="plate")

        lh.deck.assign_child_resource(tip_rack, rails=1)
        lh.deck.assign_child_resource(plate, rails=10)

        # Set initial volumes
        plate["A1"].tracker.set_liquids([(None, 200)])

        # Execute
        await lh.pick_up_tips(tip_rack["A1"])
        await lh.transfer(plate["A1"], plate["A2"], vols=100)
        await lh.drop_tips()

        # Assert
        assert plate["A1"].tracker.get_volume() == 100
        assert plate["A2"].tracker.get_volume() == 100

    finally:
        await lh.stop()
```

## 最佳实践

1. **始终首先使用仿真**：在硬件上运行之前在仿真中开发和测试协议
2. **启用跟踪**：打开尖端和体积跟踪以实现准确的可视化
3. **设置初始状态**：定义实际模拟的初始液体体积
4. **目视检查**：使用可视化工具验证平台布局和协议执行
5. **验证逻辑**：在模拟中测试边缘情况和错误条件
6. **自动化测试**：将模拟集成到 CI/CD 管道中
7. **保存布局**：使用JSON保存和共享甲板布局
8. **记录状态**：记录初始状态以实现可重复性
9. **交互式开发**：在开发过程中保持可视化工具打开
10. **协议细化**：在硬件运行之前进行仿真迭代

## 常见模式

### 开发到生产工作流程

```python
import os

# Configuration
USE_HARDWARE = os.getenv("USE_HARDWARE", "false").lower() == "true"

# Create appropriate backend
if USE_HARDWARE:
    from pylabrobot.liquid_handling.backends import STAR
    backend = STAR()
    print("Running on Hamilton STAR hardware")
else:
    from pylabrobot.liquid_handling.backends.simulation import ChatterboxBackend
    backend = ChatterboxBackend()
    print("Running in simulation mode")

# Rest of protocol is identical
lh = LiquidHandler(backend=backend, deck=STARLetDeck())

if not USE_HARDWARE:
    # Enable visualizer for simulation
    vis = Visualizer()
    await vis.start()
    lh.visualizer = vis

await lh.setup()

# Protocol execution
# ... (same code for hardware and simulation)

# Run with: USE_HARDWARE=false python protocol.py  # Simulation
# Run with: USE_HARDWARE=true python protocol.py   # Hardware
```

### 可视化协议验证

```python
async def visual_verification():
    """Run protocol with visual verification pauses"""

    vis = Visualizer()
    await vis.start()

    lh = LiquidHandler(
        backend=ChatterboxBackend(),
        deck=STARLetDeck()
    )
    lh.visualizer = vis
    await lh.setup()

    try:
        # Step 1
        await lh.pick_up_tips(tip_rack["A1:H1"])
        input("Press Enter to continue...")

        # Step 2
        await lh.aspirate(source["A1:H1"], vols=100)
        input("Press Enter to continue...")

        # Step 3
        await lh.dispense(dest["A1:H1"], vols=100)
        input("Press Enter to continue...")

        # Step 4
        await lh.drop_tips()
        input("Press Enter to finish...")

    finally:
        await lh.stop()
        await vis.stop()
```

## 故障排除

### 可视化工具未更新

- 确保在操作之前设置`lh.visualizer = vis`
- 检查是否全局启用跟踪
- 验证可视化工具正在运行 (`vis.start()`)
- 如果连接丢失，请刷新浏览器

### 跟踪不起作用

```python
# Must enable tracking BEFORE creating resources
set_tip_tracking(True)
set_volume_tracking(True)

# Then create resources
tip_rack = TIP_CAR_480_A00(name="tips")
plate = Cos_96_DW_1mL(name="plate")
```

### 模拟错误

- 模拟验证操作（例如，无法从空井中吸液）
- 使用try/except处理验证错误
- 检查初始状态是否设置正确
- 验证体积不超过容量

## 其他资源

- 可视化工具文档：https://docs.pylabrobot.org/user_guide/using-the-visualizer.html（如果有）
- 模拟指南：https://docs.pylabrobot.org/user_guide/simulation.html（如果有）
- API 参考：https://docs.pylabrobot.org/api/pylabrobot.visualizer.html
- GitHub 示例：https://github.com/PyLabRobot/pylabrobot/tree/main/examples