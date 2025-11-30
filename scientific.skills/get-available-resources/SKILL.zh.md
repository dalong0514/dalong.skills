<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：获取可用资源
描述：此技能应在任何计算密集型科学任务开始时使用，以检测和报告可用的系统资源（CPU 内核、GPU、内存、磁盘空间）。它创建一个包含资源信息和战略建议的 JSON 文件，为计算方法决策提供信息，例如是否使用并行处理（joblib、多处理）、核外计算（Dask、Zarr）、GPU 加速（PyTorch、JAX）或内存高效策略。在运行分析、训练模型、处理大型数据集或任何资源限制很重要的任务之前使用此技能。
---

# 获取可用资源

## 概述

检测可用的计算资源并为科学计算任务生成战略建议。该技能可自动识别 CPU 功能、GPU 可用性（NVIDIA CUDA、AMD ROCm、Apple Silicon Metal）、内存限制和磁盘空间，以帮助做出有关计算方法的明智决策。

## 何时使用此技能

在执行任何计算密集型任务之前主动使用此技能：

- **数据分析之前**：确定数据集是否可以加载到内存中或需要核外处理
- **模型训练之前**：检查GPU加速是否可用以及使用哪个后端
- **并行处理之前**：确定 joblib、多处理或 Dask 的最佳工作人员数量
- **大文件操作之前**：验证足够的磁盘空间和适当的存储策略
- **在项目初始化时**：了解做出架构决策的基线功能

**场景示例：**
- “帮我分析这个 50GB 的基因组数据集” → 首先使用此技能来确定是否需要 Dask/Zarr
- “根据这些数据训练神经网络” → 使用此技能来检测可用的 GPU 和后端
- “并行处理 10,000 个文件”→ 使用此技能确定最佳工作人员数量
- “运行计算密集型模拟” → 使用此技能来了解资源限制

## 此技能如何发挥作用

### 资源检测

该技能运行 `scripts/detect_resources.py` 来自动检测：

1. **CPU信息**
   - 物理和逻辑核心数
   - 处理器架构和模型
   - CPU频率信息

2. **GPU信息**
   - NVIDIA GPU：通过 nvidia-smi 检测，报告 VRAM、驱动程序版本、计算能力
   - AMD GPU：通过 rocm-smi 检测
   - Apple Silicon：检测具有金属支持和统一内存的 M1/M2/M3/M4 芯片

3. **内存信息**
   - 总内存和可用内存
   - 当前内存使用百分比
   - 交换空间可用性

4. **磁盘空间信息**
   - 工作目录的总磁盘空间和可用磁盘空间
   - 当前使用百分比

5. **操作系统信息**
   - 操作系统类型（macOS、Linux、Windows）
   - 操作系统版本和发布
   - Python版本

### 输出格式

该技能会在当前工作目录中生成一个 `.claude_resources.json` 文件，其中包含：

```json
{
  "timestamp": "2025-10-23T10:30:00",
  "os": {
    "system": "Darwin",
    "release": "25.0.0",
    "machine": "arm64"
  },
  "cpu": {
    "physical_cores": 8,
    "logical_cores": 8,
    "architecture": "arm64"
  },
  "memory": {
    "total_gb": 16.0,
    "available_gb": 8.5,
    "percent_used": 46.9
  },
  "disk": {
    "total_gb": 500.0,
    "available_gb": 200.0,
    "percent_used": 60.0
  },
  "gpu": {
    "nvidia_gpus": [],
    "amd_gpus": [],
    "apple_silicon": {
      "name": "Apple M2",
      "type": "Apple Silicon",
      "backend": "Metal",
      "unified_memory": true
    },
    "total_gpus": 1,
    "available_backends": ["Metal"]
  },
  "recommendations": {
    "parallel_processing": {
      "strategy": "high_parallelism",
      "suggested_workers": 6,
      "libraries": ["joblib", "multiprocessing", "dask"]
    },
    "memory_strategy": {
      "strategy": "moderate_memory",
      "libraries": ["dask", "zarr"],
      "note": "Consider chunking for datasets > 2GB"
    },
    "gpu_acceleration": {
      "available": true,
      "backends": ["Metal"],
      "suggested_libraries": ["pytorch-mps", "tensorflow-metal", "jax-metal"]
    },
    "large_data_handling": {
      "strategy": "disk_abundant",
      "note": "Sufficient space for large intermediate files"
    }
  }
}
```

### 战略建议

该技能会生成上下文感知的建议：

**并行处理建议：**
- **高并行性（8 个以上核心）**：使用 Dask、joblib 或多处理，workers = 核心 - 2
- **中等并行度（4-7 个核心）**：使用 joblib 或多处理，workers = cores - 1
- **顺序（< 4 核）**：首选顺序处理以避免开销

**内存策略建议：**
- **内存受限（< 4GB 可用）**：使用 Zarr、Dask 或 H5py 进行核外处理
- **中等内存（4-16GB可用）**：对于> 2GB的数据集使用Dask/Zarr
- **内存充足（> 16GB可用）**：可以将大多数数据集直接加载到内存中

**GPU 加速建议：**
- **检测到 NVIDIA GPU**：使用 PyTorch、TensorFlow、JAX、CuPy 或 RAPIDS
- **检测到 AMD GPU**：使用 PyTorch-ROCm 或 TensorFlow-ROCm
- **检测到 Apple Silicon**：将 PyTorch 与 MPS 后端、TensorFlow-Metal 或 JAX-Metal 结合使用
- **未检测到 GPU**：使用 CPU 优化的库

**大数据处理建议：**
- **磁盘受限（< 10GB）**：使用流或压缩策略
- **中等磁盘 (10-100GB)**：使用 Zarr、H5py 或 Parquet 格式
- **磁盘充足（> 100GB）**：可以自由创建大型中间文件

## 使用说明

### 第 1 步：运行资源检测

在任何计算密集型任务开始时执行检测脚本：

<<<代码块_1>>>

可选参数：
- `-o, --output <path>`：指定自定义输出路径（默认：`.claude_resources.json`）
- `-v, --verbose`：将完整资源信息打印到标准输出

### 第 2 步：阅读并应用建议

运行检测后，读取生成的 `.claude_resources.json` 文件以告知计算决策：

<<<代码块_2>>>

### 第 3 步：做出明智的决定

使用资源信息和建议做出战略选择：

**用于数据加载：**
<<<代码块_3>>>

**对于并行处理：**
<<<代码块_4>>>

**对于 GPU 加速：**
<<<代码块_5>>>

## 依赖关系

检测脚本需要以下Python包：

<<<代码块_6>>>

所有其他功能都使用 Python 标准库模块（json、os、platform、subprocess、sys、pathlib）。

## 平台支持

- **macOS**：完全支持，包括 Apple Silicon (M1/M2/M3/M4) GPU 检测
- **Linux**：完全支持包括 NVIDIA (nvidia-smi) 和 AMD (rocm-smi) GPU 检测
- **Windows**：完全支持，包括 NVIDIA GPU 检测

## 最佳实践

1. **早期运行**：在项目开始时或主要计算任务之前执行资源检测
2. **定期重新运行**：系统资源随时间变化（内存使用、磁盘空间）
3. **扩展前检查**：在扩展并行工作线程或数据大小之前验证资源
4. **记录决策**：将 `.claude_resources.json` 文件保留在项目目录中以记录资源感知决策
5. **与版本控制一起使用**：不同的机器有不同的能力；资源文件有助于保持可移植性

## 故障排除

**未检测到 GPU：**
- 确保安装 GPU 驱动程序（nvidia-smi、rocm-smi 或 Apple Silicon 的 system_profiler）
- 检查 GPU 实用程序是否在系统路径中
- 验证 GPU 没有被其他进程使用

**脚本执行失败：**
- 确保已安装 psutil：`uv pip install psutil`
- 检查Python版本兼容性（Python 3.6+）
- 验证脚本具有执行权限：`chmod +x scripts/detect_resources.py`

**内存读数不准确：**
- 内存读数是快照；实际可用内存不断变化
- 在检测准确的“可用”内存之前关闭其他应用程序
- 考虑多次运行检测并对结果取平均值