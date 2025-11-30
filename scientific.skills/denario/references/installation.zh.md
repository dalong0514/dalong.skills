<!-- 此文件由机器翻译自 installation.md -->

# 安装指南

## 系统要求

- **Python**：版本 3.12 或更高版本（必需）
- **操作系统**：Linux、macOS 或 Windows
- **虚拟环境**：建议隔离
- **LaTeX**：论文生成所需（或使用 Docker）

## 安装方法

### 方法 1：使用 uv（推荐）

uv 包管理器提供快速、可靠的依赖关系解析：

```bash
# Initialize a new project
uv init

# Add denario with app support
uv add "denario[app]"
```

### 方法 2：替代安装

使用 pip 的替代安装：

<<<代码块_1>>>

### 方法 3：从源代码构建

用于开发或定制：

<<<代码块_2>>>

### 方法四：Docker部署

Docker 提供了一个完整的环境，包含所有依赖项，包括 LaTeX：

<<<代码块_3>>>

容器启动后，通过 `http://localhost:8501` 访问 GUI。

## 验证安装

安装后，验证 denario 是否可用：

<<<代码块_4>>>

或者检查版本：

<<<代码块_5>>>

## 启动应用程序

### 命令行界面

运行图形用户界面：

<<<代码块_6>>>

这将启动一个基于 Web 的 Streamlit 应用程序，用于交互式研究工作流程管理。

### 程序化使用

直接在Python脚本中使用denario：

```python
from denario import Denario

den = Denario(project_dir="./my_project")
# Continue with workflow...
```

## 依赖关系

Denario 自动安装关键依赖项：

- **AG2**：代理编排框架
- **LangGraph**：基于图的代理工作流程
- **pandas**：数据操作
- **scikit-learn**：机器学习工具
- **matplotlib/seaborn**：可视化
- **streamlit**：GUI 框架（带有 `[app]` 额外）

## LaTeX 设置

对于纸张生成，LaTeX 必须可用：

### Linux
```bash
sudo apt-get install texlive-full
```

### macOS
```bash
brew install --cask mactex
```

### 窗口
下载并安装 [MiKTeX](https://miktex.org/download) 或 [TeX Live](https://tug.org/texlive/)。

### Docker 替代方案
Docker 镜像包含完整的 LaTeX 安装，无需手动设置。

## 安装疑难解答

### Python 版本问题

确保 Python 3.12+：
```bash
python --version
```

如果较旧，请安装较新版本或使用 pyenv 进行版本管理。

### 虚拟环境激活

**Linux/macOS：**
```bash
source venv/bin/activate
```

**Windows：**
```bash
venv\Scripts\activate
```

### 权限错误

使用 `--user` 标志或虚拟环境：
```bash
uv pip install --user "denario[app]"
```

### Docker 端口冲突

如果端口 8501 正在使用，请映射到不同的端口：
```bash
docker run -p 8502:8501 --rm pablovd/denario:latest
```

### 包冲突

创建全新的虚拟环境以避免依赖冲突。

## 更新德纳里奥

### 紫外线
```bash
uv add --upgrade denario
```

### 点
```bash
uv pip install --upgrade "denario[app]"
```

### 码头工人
```bash
docker pull pablovd/denario:latest
```

## 卸载

### 紫外线
```bash
uv remove denario
```

### 点
```bash
uv pip uninstall denario
```

### 码头工人
```bash
docker rmi pablovd/denario:latest
```