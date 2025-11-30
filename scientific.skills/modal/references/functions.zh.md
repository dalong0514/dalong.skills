<!-- 此文件由机器翻译自 functions.md -->

# 模态函数

## 基本函数定义

用 `@app.function()` 修饰 Python 函数：

```python
import modal

app = modal.App(name="my-app")

@app.function()
def my_function():
    print("Hello from Modal!")
    return "result"
```

## 调用函数

### 远程执行

调用 `.remote()` 在 Modal 上运行：

<<<代码块_1>>>

### 本地执行

调用 `.local()` 在本地运行（对于测试有用）：

<<<代码块_2>>>

## 函数参数

函数接受标准 Python 参数：

<<<代码块_3>>>

## 部署

### 临时应用程序

暂时运行：
<<<代码块_4>>>

### 已部署的应用程序

持久部署：
<<<代码块_5>>>

从其他代码访问已部署的函数：

<<<代码块_6>>>

## 入口点

### 本地入口点

在本地机器上运行的代码：

```python
@app.local_entrypoint()
def main():
    result = my_function.remote()
    print(result)
```

### 远程入口点

使用 `@app.function()` 而不使用 local_entrypoint - 完全在 Modal 上运行：

```python
@app.function()
def train_model():
    # All code runs in Modal
    ...
```

调用方式：
```bash
modal run script.py::app.train_model
```

## 参数解析

具有原始类型参数的入口点会自动进行 CLI 解析：

```python
@app.local_entrypoint()
def main(foo: int, bar: str):
    some_function.remote(foo, bar)
```

运行：
```bash
modal run script.py --foo 1 --bar "hello"
```

对于自定义解析，接受可变长度参数：

```python
import argparse

@app.function()
def train(*arglist):
    parser = argparse.ArgumentParser()
    parser.add_argument("--foo", type=int)
    args = parser.parse_args(args=arglist)
```

## 功能配置

常用参数：

```python
@app.function(
    image=my_image,           # Custom environment
    gpu="A100",               # GPU type
    cpu=2.0,                  # CPU cores
    memory=4096,              # Memory in MB
    timeout=3600,             # Timeout in seconds
    retries=3,                # Number of retries
    secrets=[my_secret],      # Environment secrets
    volumes={"/data": vol},   # Persistent storage
)
def my_function():
    ...
```

## 并行执行

### 地图

在多个输入上并行运行函数：

```python
@app.function()
def evaluate_model(x):
    return x ** 2

@app.local_entrypoint()
def main():
    inputs = list(range(100))
    for result in evaluate_model.map(inputs):
        print(result)
```

### 星图

对于具有多个参数的函数：

```python
@app.function()
def add(a, b):
    return a + b

@app.local_entrypoint()
def main():
    results = list(add.starmap([(1, 2), (3, 4)]))
    # [3, 7]
```

### 异常处理

```python
results = my_func.map(
    range(3),
    return_exceptions=True,
    wrap_returned_exceptions=False
)
# [0, 1, Exception('error')]
```

## 异步函数

定义异步函数：

```python
@app.function()
async def async_function(x: int):
    await asyncio.sleep(1)
    return x * 2

@app.local_entrypoint()
async def main():
    result = await async_function.remote.aio(42)
```

## 生成器函数

返回流结果的迭代器：

```python
@app.function()
def generate_data():
    for i in range(10):
        yield i

@app.local_entrypoint()
def main():
    for value in generate_data.remote_gen():
        print(value)
```

## 生成函数

提交后台执行的函数：

```python
@app.function()
def process_job(data):
    # Long-running job
    return result

@app.local_entrypoint()
def main():
    # Spawn without waiting
    call = process_job.spawn(data)

    # Get result later
    result = call.get(timeout=60)
```

## 程序化执行

以编程方式运行应用程序：

```python
def main():
    with modal.enable_output():
        with app.run():
            result = some_function.remote()
```

## 指定入口点

对于多个函数，指定要运行的函数：

```python
@app.function()
def f():
    print("Function f")

@app.function()
def g():
    print("Function g")
```

运行具体功能：
```bash
modal run script.py::app.f
modal run script.py::app.g
```