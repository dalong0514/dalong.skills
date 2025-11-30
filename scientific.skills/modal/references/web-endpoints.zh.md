<!-- 此文件由机器翻译自 web-endpoints.md -->

# 网络端点

## 快速入门

使用单个装饰器创建 Web 端点：

```python
image = modal.Image.debian_slim().pip_install("fastapi[standard]")

@app.function(image=image)
@modal.fastapi_endpoint()
def hello():
    return "Hello world!"
```

## 开发和部署

### 使用 `modal serve` 进行开发

<<<代码块_1>>>

创建具有实时重新加载功能的临时应用程序。端点的更改几乎立即出现。

### 使用 `modal deploy` 进行部署

<<<代码块_2>>>

创建具有稳定 URL 的持久端点。

## 简单端点

### 查询参数

<<<代码块_3>>>

致电：
<<<代码块_4>>>

### POST 请求

<<<代码块_5>>>

致电：
<<<代码块_6>>>

### Pydantic 模型

```python
from pydantic import BaseModel

class Item(BaseModel):
    name: str
    qty: int = 42

@app.function()
@modal.fastapi_endpoint(method="POST")
def process(item: Item):
    return {"processed": item.name, "quantity": item.qty}
```

## ASGI 应用程序（FastAPI、Starlette、FastHTML）

服务完整的 ASGI 应用程序：

```python
image = modal.Image.debian_slim().pip_install("fastapi[standard]")

@app.function(image=image)
@modal.concurrent(max_inputs=100)
@modal.asgi_app()
def fastapi_app():
    from fastapi import FastAPI

    web_app = FastAPI()

    @web_app.get("/")
    async def root():
        return {"message": "Hello"}

    @web_app.post("/echo")
    async def echo(request: Request):
        body = await request.json()
        return body

    return web_app
```

## WSGI 应用程序（Flask、Django）

服务同步 Web 框架：

```python
image = modal.Image.debian_slim().pip_install("flask")

@app.function(image=image)
@modal.concurrent(max_inputs=100)
@modal.wsgi_app()
def flask_app():
    from flask import Flask, request

    web_app = Flask(__name__)

    @web_app.post("/echo")
    def echo():
        return request.json

    return web_app
```

## 非 ASGI Web 服务器

对于具有自定义网络绑定的框架：

```python
@app.function()
@modal.concurrent(max_inputs=100)
@modal.web_server(8000)
def my_server():
    import subprocess
    # Must bind to 0.0.0.0, not 127.0.0.1
    subprocess.Popen("python -m http.server -d / 8000", shell=True)
```

## 流式响应

使用FastAPI的`StreamingResponse`：

```python
import time

def event_generator():
    for i in range(10):
        yield f"data: event {i}\n\n".encode()
        time.sleep(0.5)

@app.function(image=modal.Image.debian_slim().pip_install("fastapi[standard]"))
@modal.fastapi_endpoint()
def stream():
    from fastapi.responses import StreamingResponse
    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream"
    )
```

### 从模态函数流式传输

```python
@app.function(gpu="any")
def process_gpu():
    for i in range(10):
        yield f"data: result {i}\n\n".encode()
        time.sleep(1)

@app.function(image=modal.Image.debian_slim().pip_install("fastapi[standard]"))
@modal.fastapi_endpoint()
def hook():
    from fastapi.responses import StreamingResponse
    return StreamingResponse(
        process_gpu.remote_gen(),
        media_type="text/event-stream"
    )
```

### 使用 .map()

```python
@app.function()
def process_segment(i):
    return f"segment {i}\n"

@app.function(image=modal.Image.debian_slim().pip_install("fastapi[standard]"))
@modal.fastapi_endpoint()
def stream_parallel():
    from fastapi.responses import StreamingResponse
    return StreamingResponse(
        process_segment.map(range(10)),
        media_type="text/plain"
    )
```

## WebSockets

支持 `@web_server`、`@asgi_app` 和 `@wsgi_app`。每个连接维护单个函数调用。与 `@modal.concurrent` 一起使用以实现多个同时连接。

支持完整的 WebSocket 协议 (RFC 6455)。每条消息最多 2 MiB。

## 身份验证

### 代理身份验证令牌

通过 Modal 进行一流的身份验证：

```python
@app.function()
@modal.fastapi_endpoint()
def protected():
    return "authenticated!"
```

使用设置中的令牌进行保护，传入标头：
- `Modal-Key`
- `Modal-Secret`

### 不记名令牌身份验证

```python
from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials

auth_scheme = HTTPBearer()

@app.function(secrets=[modal.Secret.from_name("auth-token")])
@modal.fastapi_endpoint()
async def protected(token: HTTPAuthorizationCredentials = Depends(auth_scheme)):
    import os
    if token.credentials != os.environ["AUTH_TOKEN"]:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid token"
        )
    return "success!"
```

### 客户端IP地址

```python
from fastapi import Request

@app.function()
@modal.fastapi_endpoint()
def get_ip(request: Request):
    return f"Your IP: {request.client.host}"
```

## Web 端点 URL

### 自动生成的 URL

格式：`https://<workspace>--<app>-<function>.modal.run`

带环境后缀：`https://<workspace>-<suffix>--<app>-<function>.modal.run`

### 自定义标签

```python
@app.function()
@modal.fastapi_endpoint(label="api")
def handler():
    ...
# URL: https://workspace--api.modal.run
```

### 程序化 URL 检索

```python
@app.function()
@modal.fastapi_endpoint()
def my_endpoint():
    url = my_endpoint.get_web_url()
    return {"url": url}

# From deployed function
f = modal.Function.from_name("app-name", "my_endpoint")
url = f.get_web_url()
```

### 自定义域

适用于团队和企业计划：

```python
@app.function()
@modal.fastapi_endpoint(custom_domains=["api.example.com"])
def hello(message: str):
    return {"message": f"hello {message}"}
```

多个域：
```python
@modal.fastapi_endpoint(custom_domains=["api.example.com", "api.example.net"])
```

通配符域：
```python
@modal.fastapi_endpoint(custom_domains=["*.example.com"])
```

自动生成和更新 TLS 证书。

## 性能

### 冷启动

第一个请求可能会经历冷启动（几秒钟）。 Modal 使容器为后续请求保持活动状态。

### 缩放

- 根据流量自动缩放
- 对每个容器的多个请求使用 `@modal.concurrent`
- 超出并发限制，额外的容器会启动
- 达到最大容器数时请求队列

### 速率限制

默认值：200 个请求/秒，突发倍数为 5 秒
- 超额返回429状态码
- 联系支持人员以增加限制

### 大小限制

- 请求正文：最多 4 GiB
- 响应体：无限制
- WebSocket 消息：最多 2 MiB