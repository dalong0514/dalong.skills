<!-- 此文件由机器翻译自 monitoring.md -->

# SimPy 监控和数据收集

本指南涵盖了在 SimPy 中收集数据和监控模拟行为的技术。

## 监控策略

在实施监控之前，先定义三件事：

1. **监控什么**：进程、资源、事件或系统状态
2. **何时监控**：在变化时、每隔一定时间或在特定事件时进行监控
3. **如何存储数据**：列表、文件、数据库或实时输出

## 1. 过程监控

### 状态变量跟踪

通过记录变量变化来跟踪进程状态。

```python
import simpy

def customer(env, name, service_time, log):
    arrival_time = env.now
    log.append(('arrival', name, arrival_time))

    yield env.timeout(service_time)

    departure_time = env.now
    log.append(('departure', name, departure_time))

    wait_time = departure_time - arrival_time
    log.append(('wait_time', name, wait_time))

env = simpy.Environment()
log = []

env.process(customer(env, 'Customer 1', 5, log))
env.process(customer(env, 'Customer 2', 3, log))
env.run()

print('Simulation log:')
for entry in log:
    print(entry)
```

### 时间序列数据收集

<<<代码块_1>>>

### 多变量跟踪

<<<代码块_2>>>

## 2. 资源监控

### 猴子修补资源

修补资源方法以拦截和记录操作。

<<<代码块_3>>>

### 资源子类化

创建具有内置监控的自定义资源类。

<<<代码块_4>>>

### 容器液位监控

<<<代码块_5>>>

## 3. 事件追踪

### 环境步骤监控

通过修补环境的步骤函数来监视所有事件。

<<<代码块_6>>>

### 事件调度监视器

跟踪活动安排时间。

```python
import simpy

class MonitoredEnvironment(simpy.Environment):
    def __init__(self):
        super().__init__()
        self.scheduled_events = []

    def schedule(self, event, priority=simpy.core.NORMAL, delay=0):
        super().schedule(event, priority, delay)
        scheduled_time = self.now + delay
        self.scheduled_events.append((scheduled_time, priority, type(event).__name__))

def process(env, name, delay):
    print(f'{name}: Scheduling timeout for {delay} at {env.now}')
    yield env.timeout(delay)
    print(f'{name}: Resumed at {env.now}')

env = MonitoredEnvironment()
env.process(process(env, 'Process 1', 5))
env.process(process(env, 'Process 2', 3))
env.run()

print('\nScheduled events:')
for time, priority, event_type in env.scheduled_events:
    print(f'Time {time}, Priority {priority}, Type {event_type}')
```

## 4. 统计监控

### 队列统计

```python
import simpy

class QueueStatistics:
    def __init__(self):
        self.arrival_times = []
        self.departure_times = []
        self.queue_lengths = []
        self.wait_times = []

    def record_arrival(self, time, queue_length):
        self.arrival_times.append(time)
        self.queue_lengths.append(queue_length)

    def record_departure(self, arrival_time, departure_time):
        self.departure_times.append(departure_time)
        self.wait_times.append(departure_time - arrival_time)

    def average_wait_time(self):
        return sum(self.wait_times) / len(self.wait_times) if self.wait_times else 0

    def average_queue_length(self):
        return sum(self.queue_lengths) / len(self.queue_lengths) if self.queue_lengths else 0

def customer(env, resource, stats):
    arrival_time = env.now
    stats.record_arrival(arrival_time, len(resource.queue))

    with resource.request() as req:
        yield req
        departure_time = env.now
        stats.record_departure(arrival_time, departure_time)
        yield env.timeout(2)

env = simpy.Environment()
resource = simpy.Resource(env, capacity=1)
stats = QueueStatistics()

for i in range(5):
    env.process(customer(env, resource, stats))

env.run()

print(f'Average wait time: {stats.average_wait_time():.2f}')
print(f'Average queue length: {stats.average_queue_length():.2f}')
```

## 5. 数据导出

### CSV 导出

```python
import simpy
import csv

def export_to_csv(data, filename):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Time', 'Metric', 'Value'])
        writer.writerows(data)

def monitored_simulation(env, data_log):
    for i in range(10):
        data_log.append((env.now, 'queue_length', i % 3))
        data_log.append((env.now, 'utilization', (i % 3) / 10))
        yield env.timeout(1)

env = simpy.Environment()
data = []
env.process(monitored_simulation(env, data))
env.run()

export_to_csv(data, 'simulation_data.csv')
print('Data exported to simulation_data.csv')
```

### 实时绘图（需要 matplotlib）

```python
import simpy
import matplotlib.pyplot as plt

class RealTimePlotter:
    def __init__(self):
        self.times = []
        self.values = []

    def update(self, time, value):
        self.times.append(time)
        self.values.append(value)

    def plot(self, title='Simulation Results'):
        plt.figure(figsize=(10, 6))
        plt.plot(self.times, self.values)
        plt.xlabel('Time')
        plt.ylabel('Value')
        plt.title(title)
        plt.grid(True)
        plt.show()

def monitored_process(env, plotter):
    value = 0
    for i in range(20):
        value = value * 0.9 + (i % 5)
        plotter.update(env.now, value)
        yield env.timeout(1)

env = simpy.Environment()
plotter = RealTimePlotter()
env.process(monitored_process(env, plotter))
env.run()

plotter.plot('Process Value Over Time')
```

## 最佳实践

1. **最小化开销**：仅监控必要的内容；过多的日志记录会减慢模拟速度

2. **结构化数据**：对复杂数据点使用类或命名元组

3. **时间戳**：始终在监控数据中包含时间戳

4. **聚合**：对于长时间模拟，聚合数据而不是存储每个事件

5. **惰性评估**：考虑在模拟后收集原始数据和计算统计数据

6. **内存管理**：对于很长的模拟，定期将数据刷新到磁盘

7. **验证**：验证监控代码不会影响模拟行为

8. **关注点分离**：将监控代码与模拟逻辑分开

9. **可重用组件**：创建可在模拟中重用的通用监控类