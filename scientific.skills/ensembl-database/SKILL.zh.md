<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：ensebl-数据库
描述：“查询 250 多个物种的 Ensembl 基因组数据库 REST API。基因查找、序列检索、变异分析、比较基因组学、直向同源物、VEP 预测，用于基因组研究。”
---

# 集成数据库

## 概述

访问和查询 Ensembl 基因组数据库，这是由 EMBL-EBI 维护的脊椎动物基因组数据的综合资源。该数据库提供 250 多个物种的基因注释、序列、变异、监管信息和比较基因组学数据。当前版本为 115（2025 年 9 月）。

## 何时使用此技能

该技能应该在以下情况下使用：

- 通过符号或Ensembl ID查询基因信息
- 检索 DNA、转录本或蛋白质序列
- 使用变异效应预测器 (VEP) 分析遗传变异
- 寻找跨物种的直向同源物和旁系同源物
- 访问监管特征和基因组注释
- 基因组组装之间的坐标转换（例如，GRCh37 到 GRCh38）
- 进行比较基因组学分析
- 将 Ensembl 数据整合到基因组研究流程中

## 核心能力

### 1.基因信息检索

通过符号、Ensembl ID 或外部数据库标识符查询基因数据。

**常用操作：**
- 按符号查找基因信息（例如“BRCA2”、“TP53”）
- 检索转录本和蛋白质信息
- 获取基因坐标和染色体位置
- 访问外部数据库的交叉引用（UniProt、RefSeq 等）

**使用ensembl_rest包：**
```python
from ensembl_rest import EnsemblClient

client = EnsemblClient()

# Look up gene by symbol
gene_data = client.symbol_lookup(
    species='human',
    symbol='BRCA2'
)

# Get detailed gene information
gene_info = client.lookup_id(
    id='ENSG00000139618',  # BRCA2 Ensembl ID
    expand=True
)
```

**直接 REST API（无包）：**
<<<代码块_1>>>

### 2. 序列检索

获取各种格式（JSON、FASTA、纯文本）的基因组、转录本或蛋白质序列。

**操作：**
- 获取基因或基因组区域的 DNA 序列
- 检索转录序列（cDNA）
- 获取蛋白质序列
- 提取具有侧翼区域或修饰的序列

**示例：**
<<<代码块_2>>>

### 3. 变异分析

使用变异效应预测器 (VEP) 查询遗传变异数据并预测变异后果。

**能力：**
- 通过 rsID 或基因组坐标查找变体
- 预测变异的功能后果
- 访问人口频率数据
- 检索表型关联

**VEP 示例：**
<<<代码块_3>>>

### 4.比较基因组学

进行跨物种比较，以确定直向同源物、旁系同源物和进化关系。

**操作：**
- 查找直向同源物（不同物种中的相同基因）
- 识别旁系同源物（同一物种中的相关基因）
- 访问显示进化关系的基因树
- 检索基因家族信息

**示例：**
<<<代码块_4>>>

### 5. 基因组区域分析

查找特定区域中的所有基因组特征（基因、转录本、调控元件）。

**使用案例：**
- 识别染色体区域中的所有基因
- 寻找监管特征（启动子、增强子）
- 查找区域内的变体
- 检索结构特征

**示例：**
<<<代码块_5>>>

### 6. 装配映射

转换不同基因组组件之间的坐标（例如，GRCh37 到 GRCh38）。

**重要提示：** 使用 `https://grch37.rest.ensembl.org` 进行 GRCh37/hg19 查询，使用 `https://rest.ensembl.org` 进行当前程序集。

**示例：**
<<<代码块_6>>>

## API 最佳实践

### 速率限制

Ensembl REST API 有速率限制。请遵循以下做法：

1. **尊重率限制：** 匿名用户每秒最多 15 个请求
2. **处理 429 响应：** 当速率受限时，检查 `Retry-After` 标头并等待
3. **使用批量端点：** 查询多个项目时，在可用的情况下使用批量端点
4. **缓存结果：** 存储经常访问的数据以减少API调用

### 错误处理

始终实施正确的错误处理：

```python
import requests
import time

def query_ensembl(endpoint, params=None, max_retries=3):
    server = "https://rest.ensembl.org"
    headers = {"Content-Type": "application/json"}

    for attempt in range(max_retries):
        response = requests.get(
            f"{server}{endpoint}",
            headers=headers,
            params=params
        )

        if response.status_code == 200:
            return response.json()
        elif response.status_code == 429:
            # Rate limited - wait and retry
            retry_after = int(response.headers.get('Retry-After', 1))
            time.sleep(retry_after)
        else:
            response.raise_for_status()

    raise Exception(f"Failed after {max_retries} attempts")
```

## 安装

### Python 包（推荐）

```bash
uv pip install ensembl_rest
```

`ensembl_rest` 包为所有 Ensembl REST API 端点提供了 Pythonic 接口。

### 直接 REST API

无需安装 - 使用标准 HTTP 库，例如 `requests`：

```bash
uv pip install requests
```

## 资源

###参考资料/

- `api_endpoints.md`：所有 17 个 API 端点类别的综合文档，包含示例和参数

### 脚本/

- `ensembl_query.py`：可重用的 Python 脚本，用于常见的 Ensebl 查询，具有内置的速率限制和错误处理功能

## 常见工作流程

### 工作流程 1：基因注释管道

1. 通过符号查找基因以获得Ensembl ID
2. 检索成绩单信息
3. 获取所有转录本的蛋白质序列
4. 寻找其他物种的直向同源物
5. 导出结果

### 工作流程 2：变异分析

1.通过rsID或坐标查询变量
2. 使用 VEP 预测功能后果
3. 检查人口频率
4. 检索表型关联
5. 生成报告

### 工作流程 3：比较分析

1. 从参考物种中感兴趣的基因开始
2. 寻找目标物种的直向同源物
3. 检索所有直向同源物的序列
4. 比较基因结构和特征
5. 分析进化保守性

## 物种和组装信息

要查询可用的物种和组件：

```python
# List all available species
species_list = client.info_species()

# Get assembly information for a species
assembly_info = client.info_assembly(species='human')
```

常见物种标识符：
- 人类：`homo_sapiens` 或 `human`
- 鼠标：`mus_musculus` 或 `mouse`
- 斑马鱼：`danio_rerio` 或 `zebrafish`
- 果蝇：`drosophila_melanogaster`

## 其他资源

- **官方文档：** https://rest.ensembl.org/documentation
- **Python 包文档：** https://ensemblrest.readthedocs.io
- **EBI 培训：** https://www.ebi.ac.uk/training/online/courses/ensembl-rest-api/
- **Ensebl 浏览器：** https://useast.ensembl.org
- **GitHub 示例：** https://github.com/Ensembl/ensembl-rest/wiki