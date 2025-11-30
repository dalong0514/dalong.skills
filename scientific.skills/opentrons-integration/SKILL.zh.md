<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：opentrons-集成
描述：“用于 Flex/OT-2 机器人的实验室自动化平台。编写协议 API v2 协议、液体处理、硬件模块（加热器摇床、热循环仪）、实验室器具管理，用于自动移液工作流程。”
---

# Opentrons 集成

## 概述

Opentrons 是一个基于 Python 的实验室自动化平台，适用于 Flex 和 OT-2 机器人。编写用于液体处理、控制硬件模块（加热器摇床、热循环仪）、管理实验室器具、自动移液工作流程的协议 API v2 协议。

## 何时使用此技能

该技能应该在以下情况下使用：
- 用 Python 编写 Opentrons Protocol API v2 协议
- 在 Flex 或 OT-2 机器人上实现液体处理工作流程自动化
- 控制硬件模块（温度、磁力、加热摇床、热循环仪）
- 设置实验室器具配置和甲板布局
- 实施复杂的移液操作（连续稀释、板复制、PCR 设置）
- 管理吸头使用并优化协议效率
- 使用多通道移液器进行 96 孔板操作
- 在机器人执行之前模拟和测试协议

## 核心能力

### 1. 协议结构和元数据

每个 Opentrons 协议都遵循标准结构：

```python
from opentrons import protocol_api

# Metadata
metadata = {
    'protocolName': 'My Protocol',
    'author': 'Name <email@example.com>',
    'description': 'Protocol description',
    'apiLevel': '2.19'  # Use latest available API version
}

# Requirements (optional)
requirements = {
    'robotType': 'Flex',  # or 'OT-2'
    'apiLevel': '2.19'
}

# Run function
def run(protocol: protocol_api.ProtocolContext):
    # Protocol commands go here
    pass
```

**关键要素：**
- 从 `opentrons` 导入 `protocol_api`
- 使用协议名称、作者、描述、apiLevel 定义 `metadata` 字典
- 可选的 `requirements` 字典用于机器人类型和 API 版本
- 实现`run()`函数接收`ProtocolContext`作为参数
- 所有协议逻辑都在 `run()` 函数内部

### 2. 加载硬件

**加载仪器（移液器）：**

<<<代码块_1>>>

常见移液器名称：
- Flex：`p50_single_flex`、`p1000_single_flex`、`p50_multi_flex`、`p1000_multi_flex`
- OT-2：`p20_single_gen2`、`p300_single_gen2`、`p1000_single_gen2`、`p20_multi_gen2`、`p300_multi_gen2`

**装载实验室器具：**

<<<代码块_2>>>

**加载模块：**

<<<代码块_3>>>

### 3. 液体处理操作

**基本操作：**

<<<代码块_4>>>

**复杂操作：**

<<<代码块_5>>>

**先进技术：**

<<<代码块_6>>>

**流量控制：**

```python
# Set flow rates (µL/s)
pipette.flow_rate.aspirate = 150
pipette.flow_rate.dispense = 300
pipette.flow_rate.blow_out = 400
```

### 4. 进入井和地点

**井访问方法：**

```python
# By name
well_a1 = plate['A1']

# By index
first_well = plate.wells()[0]

# All wells
all_wells = plate.wells()  # Returns list

# By rows
rows = plate.rows()  # Returns list of lists
row_a = plate.rows()[0]  # All wells in row A

# By columns
columns = plate.columns()  # Returns list of lists
column_1 = plate.columns()[0]  # All wells in column 1

# Wells by name (dictionary)
wells_dict = plate.wells_by_name()  # {'A1': Well, 'A2': Well, ...}
```

**定位方法：**

```python
# Top of well (default: 1mm below top)
pipette.aspirate(100, well.top())
pipette.aspirate(100, well.top(z=5))  # 5mm above top

# Bottom of well (default: 1mm above bottom)
pipette.aspirate(100, well.bottom())
pipette.aspirate(100, well.bottom(z=2))  # 2mm above bottom

# Center of well
pipette.aspirate(100, well.center())
```

### 5. 硬件模块控制

**温度模块：**

```python
# Set temperature
temp_module.set_temperature(celsius=4)

# Wait for temperature
temp_module.await_temperature(celsius=4)

# Deactivate
temp_module.deactivate()

# Check status
current_temp = temp_module.temperature  # Current temperature
target_temp = temp_module.target  # Target temperature
```

**磁性模块：**

```python
# Engage (raise magnets)
mag_module.engage(height_from_base=10)  # mm from labware base

# Disengage (lower magnets)
mag_module.disengage()

# Check status
is_engaged = mag_module.status  # 'engaged' or 'disengaged'
```

**加热摇床模块：**

```python
# Set temperature
hs_module.set_target_temperature(celsius=37)

# Wait for temperature
hs_module.wait_for_temperature()

# Set shake speed
hs_module.set_and_wait_for_shake_speed(rpm=500)

# Close labware latch
hs_module.close_labware_latch()

# Open labware latch
hs_module.open_labware_latch()

# Deactivate heater
hs_module.deactivate_heater()

# Deactivate shaker
hs_module.deactivate_shaker()
```

**热循环仪模块：**

```python
# Open lid
tc_module.open_lid()

# Close lid
tc_module.close_lid()

# Set lid temperature
tc_module.set_lid_temperature(celsius=105)

# Set block temperature
tc_module.set_block_temperature(
    temperature=95,
    hold_time_seconds=30,
    hold_time_minutes=0.5,
    block_max_volume=50  # µL per well
)

# Execute profile (PCR cycling)
profile = [
    {'temperature': 95, 'hold_time_seconds': 30},
    {'temperature': 57, 'hold_time_seconds': 30},
    {'temperature': 72, 'hold_time_seconds': 60}
]
tc_module.execute_profile(
    steps=profile,
    repetitions=30,
    block_max_volume=50
)

# Deactivate
tc_module.deactivate_lid()
tc_module.deactivate_block()
```

**吸光度板读数器：**

```python
# Initialize and read
result = plate_reader.read(wavelengths=[450, 650])

# Access readings
absorbance_data = result  # Dict with wavelength keys
```

### 6. 液体追踪和标签

**定义液体：**

```python
# Define liquid types
water = protocol.define_liquid(
    name='Water',
    description='Ultrapure water',
    display_color='#0000FF'  # Hex color code
)

sample = protocol.define_liquid(
    name='Sample',
    description='Cell lysate sample',
    display_color='#FF0000'
)
```

**将液体装入孔中：**

```python
# Load liquid into specific wells
reservoir['A1'].load_liquid(liquid=water, volume=50000)  # µL
plate['A1'].load_liquid(liquid=sample, volume=100)

# Mark wells as empty
plate['B1'].load_empty()
```

### 7. 协议控制和实用程序

**执行控制：**

```python
# Pause protocol
protocol.pause(msg='Replace tip box and resume')

# Delay
protocol.delay(seconds=60)
protocol.delay(minutes=5)

# Comment (appears in logs)
protocol.comment('Starting serial dilution')

# Home robot
protocol.home()
```

**条件逻辑：**

```python
# Check if simulating
if protocol.is_simulating():
    protocol.comment('Running in simulation mode')
else:
    protocol.comment('Running on actual robot')
```

**导轨灯（仅限 Flex）：**

```python
# Turn lights on
protocol.set_rail_lights(on=True)

# Turn lights off
protocol.set_rail_lights(on=False)
```

### 8. 多通道和 8 通道移液

使用多通道移液器时：

```python
# Load 8-channel pipette
multi_pipette = protocol.load_instrument(
    'p300_multi_gen2',
    'left',
    tip_racks=[tips]
)

# Access entire column with single well reference
multi_pipette.transfer(
    volume=100,
    source=source_plate['A1'],  # Accesses entire column 1
    dest=dest_plate['A1']       # Dispenses to entire column 1
)

# Use rows() for row-wise operations
for row in plate.rows():
    multi_pipette.transfer(100, reservoir['A1'], row[0])
```

### 9. 常见协议模式

**连续稀释：**

```python
def run(protocol: protocol_api.ProtocolContext):
    # Load labware
    tips = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'D1')
    reservoir = protocol.load_labware('nest_12_reservoir_15ml', 'D2')
    plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 'D3')

    # Load pipette
    p300 = protocol.load_instrument('p300_single_flex', 'left', tip_racks=[tips])

    # Add diluent to all wells except first
    p300.transfer(100, reservoir['A1'], plate.rows()[0][1:])

    # Serial dilution across row
    p300.transfer(
        100,
        plate.rows()[0][:11],  # Source: wells 0-10
        plate.rows()[0][1:],   # Dest: wells 1-11
        mix_after=(3, 50),     # Mix 3x with 50µL after dispense
        new_tip='always'
    )
```

**板复制：**

```python
def run(protocol: protocol_api.ProtocolContext):
    # Load labware
    tips = protocol.load_labware('opentrons_flex_96_tiprack_1000ul', 'C1')
    source = protocol.load_labware('corning_96_wellplate_360ul_flat', 'D1')
    dest = protocol.load_labware('corning_96_wellplate_360ul_flat', 'D2')

    # Load pipette
    p1000 = protocol.load_instrument('p1000_single_flex', 'left', tip_racks=[tips])

    # Transfer from all wells in source to dest
    p1000.transfer(
        100,
        source.wells(),
        dest.wells(),
        new_tip='always'
    )
```

**PCR 设置：**

```python
def run(protocol: protocol_api.ProtocolContext):
    # Load thermocycler
    tc_mod = protocol.load_module('thermocyclerModuleV2')
    tc_plate = tc_mod.load_labware('nest_96_wellplate_100ul_pcr_full_skirt')

    # Load tips and reagents
    tips = protocol.load_labware('opentrons_flex_96_tiprack_200ul', 'C1')
    reagents = protocol.load_labware('opentrons_24_tuberack_nest_1.5ml_snapcap', 'D1')

    # Load pipette
    p300 = protocol.load_instrument('p300_single_flex', 'left', tip_racks=[tips])

    # Open thermocycler lid
    tc_mod.open_lid()

    # Distribute master mix
    p300.distribute(
        20,
        reagents['A1'],
        tc_plate.wells(),
        new_tip='once'
    )

    # Add samples (example for first 8 wells)
    for i, well in enumerate(tc_plate.wells()[:8]):
        p300.transfer(5, reagents.wells()[i+1], well, new_tip='always')

    # Run PCR
    tc_mod.close_lid()
    tc_mod.set_lid_temperature(105)

    # PCR profile
    tc_mod.set_block_temperature(95, hold_time_seconds=180)

    profile = [
        {'temperature': 95, 'hold_time_seconds': 15},
        {'temperature': 60, 'hold_time_seconds': 30},
        {'temperature': 72, 'hold_time_seconds': 30}
    ]
    tc_mod.execute_profile(steps=profile, repetitions=35, block_max_volume=25)

    tc_mod.set_block_temperature(72, hold_time_minutes=5)
    tc_mod.set_block_temperature(4)

    tc_mod.deactivate_lid()
    tc_mod.open_lid()
```

## 最佳实践

1. **始终指定API级别**：在元数据中使用最新的稳定API版本
2. **使用有意义的标签**：为实验室器具贴上标签，以便在日志中更容易识别
3. **检查吸头可用性**：确保有足够的吸头来完成实验方案
4. **添加注释**：使用 `protocol.comment()` 进行调试和记录
5. **首先模拟**：在机器人上运行之前始终在模拟中测试协议
6. **优雅地处理错误**：在需要时添加暂停以进行手动干预
7. **考虑时间**：当方案需要潜伏期时使用延迟
8. **跟踪液体**：使用液体跟踪来更好地验证设置
9. **优化提示使用**：在适当的时候使用 `new_tip='once'` 来保存提示
10. **控制流量**：调整粘性或挥发性液体的流量

## 故障排除

**常见问题：**

- **吸头用完**：验证吸头架容量是否符合方案要求
- **实验室器具碰撞**：检查甲板布局是否存在空间冲突
- **体积错误**：确保体积不超过孔或移液器容量
- **模块没有响应**：验证模块是否正确连接并且固件是否已更新
- **体积不准确**：校准移液器并检查是否有气泡
- **协议在模拟中失败**：检查 API 版本兼容性和实验室软件定义

## 资源

有关详细的 API 文档，请参阅此技能目录中的 `references/api_reference.md`。

有关协议模板的示例，请参阅 `scripts/` 目录。