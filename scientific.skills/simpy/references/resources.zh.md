<!-- 此文件由机器翻译自 resources.md -->

# SimPy 共享资源

本指南涵盖了 SimPy 中用于建模拥塞点和资源分配的所有资源类型。

## 资源类型概述

SimPy 提供了三大类共享资源：

1. **资源** - 有限的容量资源（例如气泵、服务器）
2. **容器** - 均质散装材料（例如燃料箱、筒仓）
3. **Stores** - Python 对象存储（例如，项目队列、仓库）

## 1. 资源

对一次可由有限数量的进程使用的模型资源。

### 资源（基本）

基本资源是具有指定容量的信号量。

```python
import simpy

env = simpy.Environment()
resource = simpy.Resource(env, capacity=2)

def process(env, resource, name):
    with resource.request() as req:
        yield req
        print(f'{name} has the resource at {env.now}')
        yield env.timeout(5)
        print(f'{name} releases the resource at {env.now}')

env.process(process(env, resource, 'Process 1'))
env.process(process(env, resource, 'Process 2'))
env.process(process(env, resource, 'Process 3'))
env.run()
```

**关键属性：**
- `capacity` - 最大并发用户数（默认值：1）
- `count` - 当前用户数
- `queue` - 排队请求列表

### 优先资源

使用优先级扩展基本资源（数字越小=优先级越高）。

<<<代码块_1>>>

**使用案例：**
- 紧急服务（救护车优先于普通车辆）
- VIP客户队列
- 具有优先级的作业调度

### 抢占资源

允许高优先级请求中断低优先级用户。

<<<代码块_2>>>

**使用案例：**
- 操作系统CPU调度
- 急诊室分诊
- 网络数据包优先级

## 2. 容器

模拟均质散装材料（连续或离散）的生产和消耗。

<<<代码块_3>>>

**关键属性：**
- `capacity` - 最大金额（默认值：float('inf')）
- `level` - 当前金额
- `init` - 初始金额（默认值：0）

**操作：**
- `put(amount)` - 添加到容器（如果已满则阻塞）
- `get(amount)` - 从容器中删除（如果不足则阻塞）

**使用案例：**
- 加油站油箱
- 制造中的缓冲存储
- 水库
- 电池电量

## 3. 商店

对 Python 对象的生产和消费进行建模。

### 商店（基本）

通用 FIFO 对象存储。

<<<代码块_4>>>

**关键属性：**
- `capacity` - 最大项目数（默认值：float('inf')）
- `items` - 存储项目列表

**操作：**
- `put(item)` - 添加项目到存储（如果已满则阻塞）
- `get()` - 删除并返回项目（如果为空则阻塞）

### 过滤器存储

允许根据过滤功能检索特定对象。

<<<代码块_5>>>

**使用案例：**
- 仓库物品拣选（特定SKU）
- 具有技能匹配的工作队列
- 按目的地进行数据包路由

### 优先存储

按优先级顺序检索的项目（最低的优先）。

<<<代码块_6>>>

**使用案例：**
- 任务调度
- 打印作业队列
- 消息优先级

## 选择正确的资源类型

|场景|资源类型 |
|----------|--------------|
|有限的服务器/机器|资源 |
|基于优先级的排队 |优先资源 |
|抢占式调度 |抢占资源 |
|燃料、水、散装材料|集装箱|
|通用项目队列 (FIFO) |商店 |
|选择性项目检索 |过滤器商店|
|优先订购商品 |优先存储 |

## 最佳实践

1. **容量规划**：根据系统约束设置实际容量
2. **请求模式**：使用上下文管理器（`with resource.request()`）进行自动清理
3. **错误处理**：将抢占式资源包装在try- except中进行中断处理
4. **监控**：跟踪队列长度和利用率（参见monitoring.md）
5. **性能**：FilterStore和PriorityStore的检索时间为O(n)；明智地用于大型商店