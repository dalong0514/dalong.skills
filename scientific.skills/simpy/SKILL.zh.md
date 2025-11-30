<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：简单
描述：Python 中基于过程的离散事件模拟框架。在构建具有流程、队列、资源和基于时间的事件的系统（例如制造系统、服务运营、网络流量、物流或实体随时间与共享资源交互的任何系统）的模拟时，请使用此技能。
---

# SimPy - 离散事件模拟

## 概述

SimPy 是一个基于标准 Python 的基于流程的离散事件仿真框架。使用 SimPy 对实体（客户、车辆、数据包等）相互交互并随着时间的推移争夺共享资源（服务器、机器、带宽等）的系统进行建模。

**核心能力：**
- 使用 Python 生成器函数进行流程建模
- 共享资源管理（服务器、容器、存储）
- 事件驱动的调度和同步
- 与挂钟时间同步的实时模拟
- 全面的监控和数据收集

## 何时使用此技能

在以下情况下使用 SimPy 技能：

1. **离散事件系统建模** - 事件以不规则间隔发生的系统
2. **资源争用** - 实体争夺有限的资源（服务器、机器、员工）
3. **队列分析** - 研究等待线、服务时间和吞吐量
4. **流程优化** - 分析制造、物流或服务流程
5. **网络模拟** - 数据包路由、带宽分配、延迟分析
6. **容量规划** - 确定实现所需性能的最佳资源水平
7. **系统验证** - 在实施之前测试系统行为

**不适合：**
- 具有固定时间步长的连续模拟（考虑 SciPy ODE 求解器）
- 独立进程，无资源共享
- 纯数学优化（考虑SciPy优化）

## 快速入门

### 基本模拟结构

```python
import simpy

def process(env, name):
    """A simple process that waits and prints."""
    print(f'{name} starting at {env.now}')
    yield env.timeout(5)
    print(f'{name} finishing at {env.now}')

# Create environment
env = simpy.Environment()

# Start processes
env.process(process(env, 'Process 1'))
env.process(process(env, 'Process 2'))

# Run simulation
env.run(until=10)
```

### 资源使用模式

<<<代码块_1>>>

## 核心概念

### 1.环境

模拟环境管理时间并安排事件。

<<<代码块_2>>>

### 2. 流程

进程是使用 Python 生成器函数（带有 `yield` 语句的函数）定义的。

<<<代码块_3>>>

### 3. 活动

事件是进程同步的基本机制。流程产生事件并在触发这些事件时恢复。

**常见事件类型：**
- `env.timeout(delay)` - 等待时间过去
- `resource.request()` - 请求资源
- `env.event()` - 创建自定义事件
- `env.process(func())` - 作为事件处理
- `event1 & event2` - 等待所有事件 (AllOf)
- `event1 | event2` - 等待任何事件 (AnyOf)

## 资源

SimPy针对不同场景提供了多种资源类型。有关全面的详细信息，请参阅`references/resources.md`。

### 资源类型摘要

|资源类型 |使用案例|
|----------------|----------|
|资源 |容量有限（服务器、机器）|
|优先资源 |基于优先级的排队 |
|抢占资源|高优先级可以打断低优先级|
|集装箱|散装材料（燃料、水）|
|商店 | Python 对象存储（先进先出）|
|过滤器商店 |选择性项目检索 |
|优先存储 |优先订购商品 |

### 快速参考

<<<代码块_4>>>

## 常见模拟模式

### 模式 1：客户-服务器队列

<<<代码块_5>>>

### 模式 2：生产者-消费者

<<<代码块_6>>>

### 模式 3：并行任务执行

```python
import simpy

def task(env, name, duration):
    print(f'{name} starting at {env.now}')
    yield env.timeout(duration)
    print(f'{name} done at {env.now}')
    return f'{name} result'

def coordinator(env):
    # Start tasks in parallel
    task1 = env.process(task(env, 'Task 1', 5))
    task2 = env.process(task(env, 'Task 2', 3))
    task3 = env.process(task(env, 'Task 3', 4))

    # Wait for all to complete
    results = yield task1 & task2 & task3
    print(f'All done at {env.now}')

env = simpy.Environment()
env.process(coordinator(env))
env.run()
```

## 工作流程指南

### 第 1 步：定义系统

识别：
- **实体**：什么在系统中移动？ （客户、零件、包裹）
- **资源**：有哪些限制？ （服务器、机器、带宽）
- **流程**：有哪些活动？ （到达、服务、出发）
- **指标**：衡量什么？ （等待时间、利用率、吞吐量）

### 第 2 步：实现流程功能

为每个流程类型创建生成器函数：

```python
def entity_process(env, name, resources, parameters):
    # Arrival logic
    arrival_time = env.now

    # Request resources
    with resource.request() as req:
        yield req

        # Service logic
        service_time = calculate_service_time(parameters)
        yield env.timeout(service_time)

    # Departure logic
    collect_statistics(env.now - arrival_time)
```

### 第 3 步：设置监控

使用监控实用程序来收集数据。有关综合技术，请参阅`references/monitoring.md`。

```python
from scripts.resource_monitor import ResourceMonitor

# Create and monitor resource
resource = simpy.Resource(env, capacity=2)
monitor = ResourceMonitor(env, resource, "Server")

# After simulation
monitor.report()
```

### 第 4 步：运行并分析

```python
# Run simulation
env.run(until=simulation_time)

# Generate reports
monitor.report()
stats.report()

# Export data for further analysis
monitor.export_csv('results.csv')
```

## 高级功能

### 进程交互

进程可以通过事件、进程产量和中断进行交互。有关详细模式，请参阅 `references/process-interaction.md`。

**关键机制：**
- **事件信号**：用于协调的共享事件
- **进程产量**：等待其他进程完成
- **中断**：强制恢复进程以抢占
### 实时模拟

将仿真与硬件在环或交互式应用程序的挂钟时间同步。请参阅`references/real-time.md`。

```python
import simpy.rt

env = simpy.rt.RealtimeEnvironment(factor=1.0)  # 1:1 time mapping
# factor=0.5 means 1 sim unit = 0.5 seconds (2x faster)
```

### 全面监控

监视进程、资源和事件。请参阅 `references/monitoring.md` 了解技术，包括：
- 状态变量跟踪
- 资源猴子修补
- 事件追踪
- 统计收集

## 脚本和模板

### basic_simulation_template.py

用于构建队列模拟的完整模板：
- 可配置参数
- 统计收集
- 客户生成
- 资源使用情况
- 报告生成

**用途：**
```python
from scripts.basic_simulation_template import SimulationConfig, run_simulation

config = SimulationConfig()
config.num_resources = 2
config.sim_time = 100
stats = run_simulation(config)
stats.report()
```

### 资源监视器.py

可重复使用的监控实用程序：
- `ResourceMonitor` - 跟踪单个资源
- `MultiResourceMonitor` - 监控多个资源
- `ContainerMonitor` - 跟踪容器级别
- 自动统计计算
- CSV 导出功能

**用途：**
```python
from scripts.resource_monitor import ResourceMonitor

monitor = ResourceMonitor(env, resource, "My Resource")
# ... run simulation ...
monitor.report()
monitor.export_csv('data.csv')
```

## 参考文档

特定主题的详细指南：

- **`references/resources.md`** - 所有资源类型及示例
- **`references/events.md`** - 事件系统和模式
- **`references/process-interaction.md`** - 进程同步
- **`references/monitoring.md`** - 数据收集技术
- **`references/real-time.md`** - 实时模拟设置

## 最佳实践

1. **生成器函数**：在处理函数中始终使用`yield`
2. **资源上下文管理器**：使用`with resource.request() as req:`进行自动清理
3. **再现性**：设置 `random.seed()` 以获得一致的结果
4. **监控**：在整个模拟过程中收集数据，而不仅仅是在结束时收集数据
5. **验证**：将简单案例与解析解进行比较
6. **文档**：注释流程逻辑和参数选择
7. **模块化设计**：独立的流程逻辑、统计和配置

## 常见陷阱

1. **忘记yield**：进程必须产生事件才能暂停
2. **事件复用**：事件只能触发一次
3. **资源泄漏**：使用上下文管理器或确保释放
4. **阻塞操作**：避免进程中Python阻塞调用
5. **时间单位**：与时间单位解释保持一致
6. **死锁**：确保至少一个进程可以取得进展

## 示例用例

- **制造**：机器调度、生产线、库存管理
- **医疗保健**：急诊室模拟、患者流程、人员分配
- **电信**：网络流量、数据包路由、带宽分配
- **交通**：交通流量、物流、车辆路线
- **服务运营**：呼叫中心、零售结账、预约安排
- **计算机系统**：CPU调度、内存管理、I/O操作