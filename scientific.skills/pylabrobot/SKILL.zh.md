<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pylabrobot
描述：实验室自动化工具包，用于控制液体处理器、读板器、泵、加热摇床、培养箱、离心机和分析设备。在自动化实验室工作流程、对液体处理机器人（Hamilton STAR、Opentrons OT-2、Tecan EVO）进行编程、集成实验室设备、管理平台布局和资源（板、吸头、容器）、读取板或创建可重复的实验室协议时，请使用此技能。适用于模拟协议和物理硬件控制。
---

#PyLab机器人

## 概述

PyLabRobot 是一款与硬件无关的纯 Python 软件开发套件，适用于自动化和自主实验室。使用此技能通过跨平台（Windows、macOS、Linux）工作的统一 Python 界面来控制液体处理机器人、读板器、泵、加热摇床、培养箱、离心机和其他实验室自动化设备。

## 何时使用此技能

在以下情况下使用此技能：
- 液体处理机器人编程（Hamilton STAR/STARlet、Opentrons OT-2、Tecan EVO）
- 自动化实验室工作流程，包括移液、样品制备或分析测量
- 管理甲板布局和实验室资源（板、吸头、容器、槽）
- 集成多个实验室设备（液体处理器、读板器、加热摇床、泵）
- 通过状态管理创建可重复的实验室协议
- 在物理硬件上运行之前模拟协议
- 使用 BMG CLARIOstar 或其他支持的读板器读取读板
- 控制温度、摇动、离心或其他材料处理操作
- 使用 Python 进行实验室自动化

## 核心能力

PyLabRobot 通过六个主要功能领域提供全面的实验室自动化，每个功能领域的详细信息都在引用/目录中：

### 1. 液体处理 (`references/liquid-handling.md`)

控制液体处理机器人来吸取、分配和转移液体。关键操作包括：
- **基本操作**：吸取、分配、在孔之间转移液体
- **吸头管理**：自动拾取、放下和跟踪移液器吸头
- **先进技术**：多通道移液、连续稀释、板复制
- **体积跟踪**：自动跟踪井中的液体体积
- **硬件支持**：Hamilton STAR/STARlet、Opentrons OT-2、Tecan EVO 等

### 2. 资源管理 (`references/resources.md`)

在分层系统中管理实验室资源：
- **资源类型**：板、吸头架、槽、管、载体和定制实验室器具
- **甲板布局**：使用坐标系将资源分配给甲板位置
- **状态管理**：跟踪尖端存在、液体量和资源状态
- **序列化**：从 JSON 文件保存和加载牌组布局和状态
- **资源发现**：通过直观的 API 访问井、吸头和容器

### 3. 硬件后端 (`references/hardware-backends.md`)

通过后端抽象连接到不同的实验室设备：
- **液体处理器**：Hamilton STAR（完全支持）、Opentrons OT-2、Tecan EVO
- **模拟**：ChatterboxBackend 用于无需硬件的协议测试
- **平台支持**：适用于 Windows、macOS、Linux 和 Raspberry Pi
- **后端切换**：通过交换后端来更改机器人，无需重写协议

### 4. 分析设备 (`references/analytical-equipment.md`)

集成读板机和分析仪器：
- **读板器**：BMG CLARIOstar 用于吸光度、发光、荧光
- **秤**：用于质量测量的梅特勒-托利多集成
- **集成模式**：将液体处理器与分析设备相结合
- **自动化工作流程**：在设备之间自动移动印版

### 5. 物料搬运 (`references/material-handling.md`)

控制环境和物料搬运设备：
- **加热摇床**：Hamilton HeaterShaker、Inheco ThermoShake
- **培养箱**：带温度控制功能的 Inheco 和 Thermo Fisher 培养箱
- **离心机**：带有铲斗定位和旋转控制的 Agilent VSpin
- **泵**：用于流体泵送操作的 Cole Parmer Masterflex
- **温度控制**：在协议期间设置和监控温度

### 6. 可视化与模拟 (`references/visualization.md`)

可视化和模拟实验室协议：
- **浏览器可视化工具**：甲板状态的实时 3D 可视化
- **模拟模式**：无需物理硬件即可测试协议
- **状态跟踪**：视觉监控吸头存在和液体体积
- **甲板编辑器**：用于设计甲板布局的图形工具
- **协议验证**：在硬件上运行之前验证协议

## 快速入门

要开始使用 PyLabRobot，请安装软件包并初始化液体处理程序：

```python
# Install PyLabRobot
# uv pip install pylabrobot

# Basic liquid handling setup
from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends import STAR
from pylabrobot.resources import STARLetDeck

# Initialize liquid handler
lh = LiquidHandler(backend=STAR(), deck=STARLetDeck())
await lh.setup()

# Basic operations
await lh.pick_up_tips(tip_rack["A1:H1"])
await lh.aspirate(plate["A1"], vols=100)
await lh.dispense(plate["A2"], vols=100)
await lh.drop_tips()
```

## 使用参考

此技能可以组织多个参考文件中的详细信息。在以下情况下加载相关参考：
- **液体处理**：编写移液协议、吸头管理、转移
- **资源**：定义平台布局、管理板/吸头、定制实验室器具
- **硬件后端**：连接到特定机器人，切换平台
- **分析设备**：集成读板机、天平或分析设备
- **材料处理**：使用加热振荡器、培养箱、离心机、泵
- **可视化**：模拟协议，可视化甲板状态

所有参考文件都可以在 `references/` 目录中找到，并包含全面的示例、API 使用模式和最佳实践。

## 最佳实践

使用 PyLabRobot 创建实验室自动化协议时：

1. **从模拟开始**：在硬件上运行之前使用 ChatterboxBackend 和可视化工具测试协议
2. **启用跟踪**：打开笔尖跟踪和音量跟踪以进行准确的状态管理
3. **资源命名**：为所有资源（板、吸头架、容器）使用清晰的描述性名称
4. **状态序列化**：将牌组布局和状态保存为 JSON 以实现可重复性
5. **错误处理**：为硬件操作实现适当的异步错误处理
6. **温度控制**：尽早设定温度，因为加热/冷却需要时间
7. **模块化协议**：将复杂的工作流程分解为可重用的功能
8. **文档**：参考官方文档 https://docs.pylabrobot.org 了解最新功能

## 常见工作流程

### 液体传输协议

<<<代码块_1>>>

### 读板工作流程

<<<代码块_2>>>

## 其他资源

- **官方文档**：https://docs.pylabrobot.org
- **GitHub 存储库**：https://github.com/PyLabRobot/pylabrobot
- **社区论坛**：https://discuss.pylabrobot.org
- **PyPI 包**：https://pypi.org/project/PyLabRobot/

具体功能的详细使用请参考`references/`目录下对应的参考文件。