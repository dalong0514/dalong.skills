<!-- 此文件由机器翻译自 images.md -->

# 模态图像

## 概述

模态图像定义了代码在安装了依赖项的容器中运行的环境。图像是从基础图像开始通过方法链构建的。

## 基础镜像

从基础镜像和链方法开始：

```python
image = (
    modal.Image.debian_slim(python_version="3.13")
    .apt_install("git")
    .uv_pip_install("torch<3")
    .env({"HALT_AND_CATCH_FIRE": "0"})
    .run_commands("git clone https://github.com/modal-labs/agi")
)
```

可用的基础图像：
- `Image.debian_slim()` - 带有 Python 的 Debian Linux
- `Image.micromamba()` - 基于 Micromamba 包管理器
- `Image.from_registry()` - 从 Docker Hub、ECR 等拉取。
- `Image.from_dockerfile()` - 从现有 Dockerfile 构建

## 安装 Python 包

### 带紫外线（推荐）

使用 `.uv_pip_install()` 进行快速包安装：

<<<代码块_1>>>

### 用点

如果需要，回退到标准点：

<<<代码块_2>>>

紧密地固定依赖关系（例如，`"torch==2.8.0"`）以实现可重复性。

## 安装系统包

使用 apt 安装 Linux 软件包：

<<<代码块_3>>>

## 设置环境变量

将字典传递给 `.env()`：

<<<代码块_4>>>

## 运行 Shell 命令

在镜像构建过程中执行命令：

<<<代码块_5>>>

## 在构建时运行 Python 函数

下载模型权重或执行设置：

<<<代码块_6>>>

## 添加本地文件

### 添加文件或目录

```python
image = modal.Image.debian_slim().add_local_dir(
    "/user/erikbern/.aws",
    remote_path="/root/.aws"
)
```

默认情况下，文件在容器启动时添加。使用 `copy=True` 包含在构建的图像中。

### 添加Python源

添加可导入的 Python 模块：

```python
image = modal.Image.debian_slim().add_local_python_source("local_module")

@app.function(image=image)
def f():
    import local_module
    local_module.do_stuff()
```

## 使用现有容器镜像

### 来自公共登记处

```python
sklearn_image = modal.Image.from_registry("huanjason/scikit-learn")

@app.function(image=sklearn_image)
def fit_knn():
    from sklearn.neighbors import KNeighborsClassifier
    ...
```

可以从 Docker Hub、Nvidia NGC、AWS ECR、GitHub ghcr.io 拉取。

### 来自私人登记处

使用模态机密进行身份验证：

**Docker 中心**：
```python
secret = modal.Secret.from_name("my-docker-secret")
image = modal.Image.from_registry(
    "private-repo/image:tag",
    secret=secret
)
```

**AWS ECR**：
```python
aws_secret = modal.Secret.from_name("my-aws-secret")
image = modal.Image.from_aws_ecr(
    "000000000000.dkr.ecr.us-east-1.amazonaws.com/my-private-registry:latest",
    secret=aws_secret,
)
```

### 来自 Dockerfile

```python
image = modal.Image.from_dockerfile("Dockerfile")

@app.function(image=image)
def fit():
    import sklearn
    ...
```

导入后仍可使用其他图像方法进行扩展。

## 使用 Micromamba

对于Python和系统包的协调安装：

```python
numpyro_pymc_image = (
    modal.Image.micromamba()
    .micromamba_install("pymc==5.10.4", "numpyro==0.13.2", channels=["conda-forge"])
)
```

## 构建时的 GPU 支持

在 GPU 实例上运行构建步骤：

```python
image = (
    modal.Image.debian_slim()
    .pip_install("bitsandbytes", gpu="H100")
)
```

## 图像缓存

图像按层缓存。破坏一层上的缓存会导致后续层的级联重建。

最后定义经常更改的层，以最大限度地提高缓存重用率。

### 强制重建

```python
image = (
    modal.Image.debian_slim()
    .apt_install("git")
    .pip_install("slack-sdk", force_build=True)
)
```

或者设置环境变量：
```bash
MODAL_FORCE_BUILD=1 modal run ...
```

## 处理不同的本地/远程包

导入仅在函数体内远程可用的包：

```python
@app.function(image=image)
def my_function():
    import pandas as pd  # Only imported remotely
    df = pd.DataFrame()
    ...
```

或者使用导入上下文管理器：

```python
pandas_image = modal.Image.debian_slim().pip_install("pandas")

with pandas_image.imports():
    import pandas as pd

@app.function(image=pandas_image)
def my_function():
    df = pd.DataFrame()
```

## 使用 eStargz 从注册表快速拉取

通过 eStargz 压缩提高拉取性能：

```bash
docker buildx build --tag "<registry>/<namespace>/<repo>:<version>" \
  --output type=registry,compression=estargz,force-compression=true,oci-mediatypes=true \
  .
```

支持的注册表：
- AWS ECR
- Docker 中心
- 谷歌工件注册表