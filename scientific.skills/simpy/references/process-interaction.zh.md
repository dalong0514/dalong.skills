<!-- 此文件由机器翻译自 process-interaction.md -->

# SimPy 进程交互

本指南涵盖了 SimPy 模拟中进程交互和同步的机制。

## 交互机制概述

SimPy 为进程交互提供了三种主要方式：

1. **基于事件的钝化/重新激活** - 用于信令的共享事件
2. **等待进程终止** - 产生进程对象
3. **中断** - 强制恢复暂停的进程

## 1. 基于事件的钝化和重新激活

进程可以共享事件来协调它们的执行。

### 基本信号模式

```python
import simpy

def controller(env, signal_event):
    print(f'Controller: Preparing at {env.now}')
    yield env.timeout(5)
    print(f'Controller: Sending signal at {env.now}')
    signal_event.succeed()

def worker(env, signal_event):
    print(f'Worker: Waiting for signal at {env.now}')
    yield signal_event
    print(f'Worker: Received signal, starting work at {env.now}')
    yield env.timeout(3)
    print(f'Worker: Work complete at {env.now}')

env = simpy.Environment()
signal = env.event()
env.process(controller(env, signal))
env.process(worker(env, signal))
env.run()
```

**使用案例：**
- 协调行动的启动信号
- 完成通知
- 广播状态变化

### 多个服务员

多个进程可以等待同一个信号事件。

<<<代码块_1>>>

### 屏障同步

<<<代码块_2>>>

## 2.等待进程终止

进程本身就是事件，因此您可以让它们等待完成。

### 顺序流程执行

<<<代码块_3>>>

### 并行进程执行

<<<代码块_4>>>

### 第一个完成模式

<<<代码块_5>>>

## 3. 进程中断

可以使用 `process.interrupt()` 中断进程，这会引发 `Interrupt` 异常。

### 基本中断

<<<代码块_6>>>

### 可恢复中断

进程可以在中断后重新产生相同的事件以继续等待。

```python
import simpy

def resumable_worker(env):
    work_left = 10

    while work_left > 0:
        try:
            print(f'Worker: Working ({work_left} units left) at {env.now}')
            start = env.now
            yield env.timeout(work_left)
            work_left = 0
            print(f'Worker: Completed at {env.now}')
        except simpy.Interrupt:
            work_left -= (env.now - start)
            print(f'Worker: Interrupted! {work_left} units left at {env.now}')

def interrupter(env, worker_proc):
    yield env.timeout(3)
    worker_proc.interrupt()
    yield env.timeout(2)
    worker_proc.interrupt()

env = simpy.Environment()
worker_proc = env.process(resumable_worker(env))
env.process(interrupter(env, worker_proc))
env.run()
```

### 自定义原因中断

```python
import simpy

def machine(env, name):
    while True:
        try:
            print(f'{name}: Operating at {env.now}')
            yield env.timeout(5)
        except simpy.Interrupt as interrupt:
            if interrupt.cause == 'maintenance':
                print(f'{name}: Maintenance required at {env.now}')
                yield env.timeout(2)
                print(f'{name}: Maintenance complete at {env.now}')
            elif interrupt.cause == 'emergency':
                print(f'{name}: Emergency stop at {env.now}')
                break

def maintenance_scheduler(env, machine_proc):
    yield env.timeout(7)
    machine_proc.interrupt(cause='maintenance')
    yield env.timeout(10)
    machine_proc.interrupt(cause='emergency')

env = simpy.Environment()
machine_proc = env.process(machine(env, 'Machine 1'))
env.process(maintenance_scheduler(env, machine_proc))
env.run()
```

### 带中断的抢占资源

```python
import simpy

def user(env, name, resource, priority, duration):
    with resource.request(priority=priority) as req:
        try:
            yield req
            print(f'{name} (priority {priority}): Got resource at {env.now}')
            yield env.timeout(duration)
            print(f'{name}: Done at {env.now}')
        except simpy.Interrupt:
            print(f'{name}: Preempted at {env.now}')

env = simpy.Environment()
resource = simpy.PreemptiveResource(env, capacity=1)

env.process(user(env, 'Low priority user', resource, priority=10, duration=10))
env.process(user(env, 'High priority user', resource, priority=1, duration=5))
env.run()
```

## 高级模式

### 生产者-消费者与信令

```python
import simpy

class Buffer:
    def __init__(self, env, capacity):
        self.env = env
        self.capacity = capacity
        self.items = []
        self.item_available = env.event()

    def put(self, item):
        if len(self.items) < self.capacity:
            self.items.append(item)
            if not self.item_available.triggered:
                self.item_available.succeed()
            return True
        return False

    def get(self):
        if self.items:
            return self.items.pop(0)
        return None

def producer(env, buffer):
    item_id = 0
    while True:
        yield env.timeout(2)
        item = f'Item {item_id}'
        if buffer.put(item):
            print(f'Producer: Added {item} at {env.now}')
            item_id += 1

def consumer(env, buffer):
    while True:
        if buffer.items:
            item = buffer.get()
            print(f'Consumer: Retrieved {item} at {env.now}')
            yield env.timeout(3)
        else:
            print(f'Consumer: Waiting for items at {env.now}')
            yield buffer.item_available
            buffer.item_available = env.event()

env = simpy.Environment()
buffer = Buffer(env, capacity=5)
env.process(producer(env, buffer))
env.process(consumer(env, buffer))
env.run(until=20)
```

### 握手协议

```python
import simpy

def sender(env, request_event, acknowledge_event):
    for i in range(3):
        print(f'Sender: Sending request {i} at {env.now}')
        request_event.succeed(value=f'Request {i}')
        yield acknowledge_event
        print(f'Sender: Received acknowledgment at {env.now}')

        # Reset events for next iteration
        request_event = env.event()
        acknowledge_event = env.event()
        yield env.timeout(1)

def receiver(env, request_event, acknowledge_event):
    for i in range(3):
        request = yield request_event
        print(f'Receiver: Got {request} at {env.now}')
        yield env.timeout(2)  # Process request
        acknowledge_event.succeed()
        print(f'Receiver: Sent acknowledgment at {env.now}')

        # Reset for next iteration
        request_event = env.event()
        acknowledge_event = env.event()

env = simpy.Environment()
request = env.event()
ack = env.event()
env.process(sender(env, request, ack))
env.process(receiver(env, request, ack))
env.run()
```

## 最佳实践

1. **选择正确的机制**：
   - 使用事件进行信号和广播
   - 使用顺序/并行工作流程的流程产量
   - 使用中断进行抢占和紧急处理

2. **异常处理**：始终将容易中断的代码包装在 try- except 块中

3. **事件生命周期**：记住事件只能触发一次；为重复信号创建新事件

4. **进程引用**：存储进程对象，以便稍后需要中断它们

5. **原因信息**：使用中断原因来传达中断发生的原因

6. **可恢复模式**：跟踪进度以在中断后恢复

7. **避免死锁**：确保任何时候至少有一个进程可以取得进展