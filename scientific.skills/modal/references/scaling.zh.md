<!-- 此文件由机器翻译自 scaling.md -->

# 在 Modal 上横向扩展

## 自动缩放

每个模态函数都对应一个自动缩放容器池。 Modal 的自动缩放器：
- 当没有可用容量时旋转容器
- 当资源空闲时关闭容器
- 当没有输入要处理时默认缩放为零

自动缩放决策是快速而频繁地做出的。

## 使用 `.map()` 并行执行

使用不同的输入并行重复运行函数：

```python
@app.function()
def evaluate_model(x):
    return x ** 2

@app.local_entrypoint()
def main():
    inputs = list(range(100))
    # Runs 100 inputs in parallel across containers
    for result in evaluate_model.map(inputs):
        print(result)
```

### 带有 `.starmap()` 的多个参数

对于具有多个参数的函数：

<<<代码块_1>>>

### 异常处理

<<<代码块_2>>>

## 自动缩放配置

使用参数配置自动缩放器行为：

<<<代码块_3>>>

参数：
- **max_containers**：容器总数的上限
- **min_containers**：即使在不活动时也能保持最低温度
- **buffer_containers**：功能活动时的缓冲区大小（额外的输入不需要排队）
- **scaledown_window**：缩小之前的最大空闲时间（秒）

权衡：
- 更大的暖池/缓冲区→更高的成本，更低的延迟
- 更长的缩小窗口→减少不频繁请求的流失

## 动态自动缩放器更新

更新自动缩放器设置而不重新部署：

<<<代码块_4>>>

设置在下次部署时恢复为装饰器配置，或者被进一步更新覆盖：

<<<代码块_5>>>

### 基于时间的缩放

根据一天中的时间调整温水池：

<<<代码块_6>>>

### 对于班级

更新特定参数实例的自动缩放器：

```python
MyClass = modal.Cls.from_name("my-app", "MyClass")
obj = MyClass(model_version="3.5")
obj.update_autoscaler(buffer_containers=2)  # type: ignore
```

## 输入并发

使用 `@modal.concurrent` 处理每个容器的多个输入：

```python
@app.function()
@modal.concurrent(max_inputs=100)
def my_function(input: str):
    # Container can handle up to 100 concurrent inputs
    ...
```

I/O 密集型工作负载的理想选择：
- 数据库查询
- 外部API请求
- 远程模态函数调用

### 并发机制

**同步函数**：单独的线程（必须是线程安全的）

```python
@app.function()
@modal.concurrent(max_inputs=10)
def sync_function():
    time.sleep(1)  # Must be thread-safe
```

**异步函数**：单独的异步任务（不得阻止事件循环）

```python
@app.function()
@modal.concurrent(max_inputs=10)
async def async_function():
    await asyncio.sleep(1)  # Must not block event loop
```

### 目标与最大输入

```python
@app.function()
@modal.concurrent(
    max_inputs=120,    # Hard limit
    target_inputs=100  # Autoscaler target
)
def my_function(input: str):
    # Allow 20% burst above target
    ...
```

Autoscaler 的目标是 `target_inputs`，但容器在扩展期间可能会爆发到 `max_inputs`。

## 缩放限制

模态强制每个函数的限制：
- 2,000 个待处理输入（尚未分配给容器）
- 总共 25,000 个输入（正在运行 + 待处理）

对于 `.spawn()` 异步作业：最多 100 万个待处理输入。

超出限制将返回 `Resource Exhausted` 错误 - 稍后重试。

每个 `.map()` 调用：最多 1,000 个并发输入。

## 异步用法

使用异步 API 实现任意并行执行模式：

```python
@app.function()
async def async_task(x):
    await asyncio.sleep(1)
    return x * 2

@app.local_entrypoint()
async def main():
    tasks = [async_task.remote.aio(i) for i in range(100)]
    results = await asyncio.gather(*tasks)
```

## 常见问题

**不正确**：使用Python的内置地图（按顺序运行）
```python
# DON'T DO THIS
results = map(evaluate_model, inputs)
```

**不正确**：首先调用函数
```python
# DON'T DO THIS
results = evaluate_model(inputs).map()
```

**正确**：在 Modal 函数对象上调用 .map()
```python
# DO THIS
results = evaluate_model.map(inputs)
```