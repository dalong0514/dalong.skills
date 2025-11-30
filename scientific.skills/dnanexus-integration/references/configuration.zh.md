<!-- 此文件由机器翻译自 configuration.md -->

# DNAnexus 应用程序配置和依赖项

## 概述

本指南涵盖通过 dxapp.json 元数据配置应用程序以及管理依赖项（包括系统包、Python 库和 Docker 容器）。

## dxapp.json 结构

`dxapp.json` 文件是 DNAnexus 应用程序和小程序的配置文件。它定义元数据、输入、输出、执行要求和依赖项。

### 最小示例

```json
{
  "name": "my-app",
  "title": "My Analysis App",
  "summary": "Performs analysis on input files",
  "dxapi": "1.0.0",
  "version": "1.0.0",
  "inputSpec": [],
  "outputSpec": [],
  "runSpec": {
    "interpreter": "python3",
    "file": "src/my-app.py",
    "distribution": "Ubuntu",
    "release": "24.04"
  }
}
```

## 元数据字段

### 必填字段

<<<代码块_1>>>

### 可选元数据

<<<代码块_2>>>

## 输入规范

定义输入参数：

<<<代码块_3>>>

### 输入类

- `file` - 文件对象
- `record` - 记录对象
- `applet` - 小程序参考
- `string` - 文本字符串
- `int` - 整数
- `float` - 浮点数
- `boolean` - 真/假
- `hash` - 键值映射
- `array:class` - 指定类的数组

### 输入选项

- `name` - 参数名称（必填）
- `class` - 数据类型（必需）
- `optional` - 参数是否可选（默认：false）
- `default` - 可选参数的默认值
- `label` - 在 UI 中显示名称
- `help` - 描述文本
- `patterns` - 文件名模式（对于文件）
- `suggestions` - 预定义参考数据
- `choices` - 允许的值（对于字符串/数字）
- `group` - UI 分组

## 输出规格

定义输出参数：

<<<代码块_4>>>

## 运行规范

定义应用程序如何执行：

<<<代码块_5>>>

## 系统要求

### 实例类型选择

<<<代码块_6>>>

**常见实例类型**：
- `mem1_ssd1_v2_x2` - 2 核，3.9 GB RAM
- `mem1_ssd1_v2_x4` - 4 核，7.8 GB RAM
- `mem2_ssd1_v2_x4` - 4 核，15.6 GB RAM
- `mem2_ssd1_v2_x8` - 8 核，31.2 GB RAM
- `mem3_ssd1_v2_x8` - 8 核，62.5 GB RAM
- `mem3_ssd1_v2_x16` - 16 核，125 GB RAM

### 集群规格

对于分布式计算：

```json
{
  "systemRequirements": {
    "main": {
      "clusterSpec": {
        "type": "spark",
        "version": "3.1.2",
        "initialInstanceCount": 3,
        "instanceType": "mem1_ssd1_v2_x4",
        "bootstrapScript": "bootstrap.sh"
      }
    }
  }
}
```

## 区域选项

跨区域部署应用程序：

```json
{
  "regionalOptions": {
    "aws:us-east-1": {
      "systemRequirements": {
        "*": {"instanceType": "mem2_ssd1_v2_x4"}
      },
      "assetDepends": [
        {"id": "record-xxxx"}
      ]
    },
    "azure:westus": {
      "systemRequirements": {
        "*": {"instanceType": "azure:mem2_ssd1_x4"}
      }
    }
  }
}
```

## 依赖管理

### 系统包 (execDepends)

在运行时安装 Ubuntu 软件包：

```json
{
  "runSpec": {
    "execDepends": [
      {"name": "samtools"},
      {"name": "bwa"},
      {"name": "python3-pip"},
      {"name": "r-base", "version": "4.0.0"}
    ]
  }
}
```

使用 Ubuntu 存储库中的 `apt-get` 安装软件包。

### Python 依赖项

#### 选项 1：通过 execDepends 中的 pip 安装

```json
{
  "runSpec": {
    "execDepends": [
      {"name": "python3-pip"}
    ]
  }
}
```

然后在您的应用程序脚本中：
```python
import subprocess
subprocess.check_call(["pip", "install", "numpy==1.24.0", "pandas==2.0.0"])
```

#### 选项 2：需求文件

创建`resources/requirements.txt`：
```
numpy==1.24.0
pandas==2.0.0
scikit-learn==1.3.0
```

在您的应用程序中：
```python
subprocess.check_call(["pip", "install", "-r", "requirements.txt"])
```

### 捆绑依赖项

在应用程序中包含自定义工具或库：

**文件结构**：
```
my-app/
├── dxapp.json
├── src/
│   └── my-app.py
└── resources/
    ├── tools/
    │   └── custom_tool
    └── scripts/
        └── helper.py
```

访问应用程序中的资源：
```python
import os

# Resources are in parent directory
resources_dir = os.path.join(os.path.dirname(__file__), "..", "resources")
tool_path = os.path.join(resources_dir, "tools", "custom_tool")

# Run bundled tool
subprocess.check_call([tool_path, "arg1", "arg2"])
```

### 资产依赖

资产是预先构建的依赖项捆绑包，可以在应用程序之间共享。

#### 使用资产

```json
{
  "runSpec": {
    "assetDepends": [
      {
        "name": "bwa-asset",
        "id": {"$dnanexus_link": "record-xxxx"}
      }
    ]
  }
}
```

资产在运行时安装并通过环境变量访问：
```python
import os
asset_dir = os.environ.get("DX_ASSET_BWA")
bwa_path = os.path.join(asset_dir, "bin", "bwa")
```

#### 创建资产

创建资产目录：
```bash
mkdir bwa-asset
cd bwa-asset
# Install software
./configure --prefix=$PWD/usr/local
make && make install
```

建立资产：
```bash
dx build_asset bwa-asset --destination=project-xxxx:/assets/
```

## Docker 集成

### 使用 Docker 镜像

```json
{
  "runSpec": {
    "interpreter": "python3",
    "file": "src/my-app.py",
    "distribution": "Ubuntu",
    "release": "24.04",
    "systemRequirements": {
      "*": {
        "instanceType": "mem2_ssd1_v2_x4"
      }
    },
    "execDepends": [
      {"name": "docker.io"}
    ]
  }
}
```

在应用程序中使用 Docker：
```python
import subprocess

# Pull Docker image
subprocess.check_call(["docker", "pull", "biocontainers/samtools:v1.9"])

# Run command in container
subprocess.check_call([
    "docker", "run",
    "-v", f"{os.getcwd()}:/data",
    "biocontainers/samtools:v1.9",
    "samtools", "view", "/data/input.bam"
])
```

### Docker 作为基础镜像

对于完全在 Docker 中运行的应用程序：

```json
{
  "runSpec": {
    "interpreter": "bash",
    "file": "src/wrapper.sh",
    "distribution": "Ubuntu",
    "release": "24.04",
    "execDepends": [
      {"name": "docker.io"}
    ]
  }
}
```

## 访问要求

请求特殊权限：

```json
{
  "access": {
    "network": ["*"],           // Internet access
    "project": "CONTRIBUTE",    // Project write access
    "allProjects": "VIEW",      // Read other projects
    "developer": true           // Advanced permissions
  }
}
```

**网络访问**：
- `["*"]` - 完整的互联网
- `["github.com", "pypi.org"]` - 特定域

## 超时配置

```json
{
  "runSpec": {
    "timeoutPolicy": {
      "*": {
        "days": 1,
        "hours": 12,
        "minutes": 30
      }
    }
  }
}
```

## 示例：完成 dxapp.json

```json
{
  "name": "rna-seq-pipeline",
  "title": "RNA-Seq Analysis Pipeline",
  "summary": "Aligns RNA-seq reads and quantifies gene expression",
  "description": "Comprehensive RNA-seq pipeline using STAR aligner and featureCounts",
  "version": "1.0.0",
  "dxapi": "1.0.0",
  "categories": ["Read Mapping", "RNA-Seq"],

  "inputSpec": [
    {
      "name": "reads",
      "label": "FASTQ reads",
      "class": "array:file",
      "patterns": ["*.fastq.gz", "*.fq.gz"],
      "help": "Single-end or paired-end RNA-seq reads"
    },
    {
      "name": "reference_genome",
      "label": "Reference genome",
      "class": "file",
      "patterns": ["*.fa", "*.fasta"],
      "suggestions": [
        {
          "name": "Human GRCh38",
          "project": "project-reference",
          "path": "/genomes/GRCh38.fa"
        }
      ]
    },
    {
      "name": "gtf_file",
      "label": "Gene annotation (GTF)",
      "class": "file",
      "patterns": ["*.gtf", "*.gtf.gz"]
    }
  ],

  "outputSpec": [
    {
      "name": "aligned_bam",
      "label": "Aligned reads (BAM)",
      "class": "file",
      "patterns": ["*.bam"]
    },
    {
      "name": "counts",
      "label": "Gene counts",
      "class": "file",
      "patterns": ["*.counts.txt"]
    },
    {
      "name": "qc_report",
      "label": "QC report",
      "class": "file",
      "patterns": ["*.html"]
    }
  ],

  "runSpec": {
    "interpreter": "python3",
    "file": "src/rna-seq-pipeline.py",
    "distribution": "Ubuntu",
    "release": "24.04",

    "execDepends": [
      {"name": "python3-pip"},
      {"name": "samtools"},
      {"name": "subread"}
    ],

    "assetDepends": [
      {
        "name": "star-aligner",
        "id": {"$dnanexus_link": "record-star-asset"}
      }
    ],

    "systemRequirements": {
      "main": {
        "instanceType": "mem3_ssd1_v2_x16"
      }
    },

    "timeoutPolicy": {
      "*": {"hours": 8}
    }
  },

  "access": {
    "network": ["*"]
  },

  "details": {
    "contactEmail": "support@example.com",
    "upstreamVersion": "STAR 2.7.10a, Subread 2.0.3",
    "citations": ["doi:10.1093/bioinformatics/bts635"]
  }
}
```

## 最佳实践

1. **版本管理**：对应用程序使用语义版本控制
2. **实例类型**：从较小的实例开始，根据需要扩展
3. **依赖关系**：清楚地记录所有依赖关系
4. **错误消息**：为无效输入提供有用的错误消息
5. **测试**：使用各种输入类型和大小进行测试
6. **文档**：编写清晰的描述和帮助文本
7. **资源**：捆绑常用工具，避免重复下载
8. **Docker**：使用 Docker 处理复杂的依赖链
9. **资产**：为跨应用程序共享的重依赖项创建资产
10. **超时**：根据预期运行时间设置合理的超时
11. **网络访问**：仅请求必要的网络权限
12. **区域支持**：对多区域应用程序使用regionalOptions

## 常见模式

### 生物信息学工具
```json
{
  "inputSpec": [
    {"name": "input_file", "class": "file", "patterns": ["*.bam"]},
    {"name": "threads", "class": "int", "default": 4, "optional": true}
  ],
  "runSpec": {
    "execDepends": [{"name": "tool-name"}],
    "systemRequirements": {
      "main": {"instanceType": "mem2_ssd1_v2_x8"}
    }
  }
}
```

### Python数据分析

```json
{
  "runSpec": {
    "interpreter": "python3",
    "execDepends": [
      {"name": "python3-pip"}
    ],
    "systemRequirements": {
      "main": {"instanceType": "mem2_ssd1_v2_x4"}
    }
  }
}
```

### 基于 Docker 的应用程序

```json
{
  "runSpec": {
    "interpreter": "bash",
    "execDepends": [
      {"name": "docker.io"}
    ],
    "systemRequirements": {
      "main": {"instanceType": "mem2_ssd1_v2_x8"}
    }
  },
  "access": {
    "network": ["*"]
  }
}
```