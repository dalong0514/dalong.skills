<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：反应组数据库
描述：“查询 Reactome REST API，用于系统生物学研究的通路分析、富集、基因通路图谱、疾病通路、分子相互作用、表达分析。”
---

# 反应组数据库

## 概述

Reactome 是一个免费、开源、精心策划的通路数据库，包含 2,825 多个人类通路。查询生物通路、执行过度表达和表达分析、将基因映射到通路、通过 REST API 和 Python 客户端探索分子相互作用以进行系统生物学研究。

## 何时使用此技能

该技能应该在以下情况下使用：
- 对基因或蛋白质列表进行通路富集分析
- 分析基因表达数据以确定相关的生物学途径
- 查询特定途径信息、反应或分子相互作用
- 将基因或蛋白质映射到生物途径和过程
- 探索疾病相关途径和机制
- 在 Reactome Pathway Browser 中可视化分析结果
- 进行跨物种比较路径分析

## 核心能力

Reactome 提供两个主要的 API 服务和一个 Python 客户端库：

### 1.内容服务-数据检索

查询和检索生物途径数据、分子相互作用和实体信息。

**常用操作：**
- 检索路径信息和层次结构
- 查询特定实体（蛋白质、反应、复合物）
- 获取途径中的参与分子
- 访问数据库版本和元数据
- 探索通道隔间和位置

**API 基本 URL：** `https://reactome.org/ContentService`

### 2.分析服务-通路分析

对基因列表和表达数据进行计算分析。

**分析类型：**
- **过度代表性分析**：从基因/蛋白质列表中识别具有统计意义的途径
- **表达数据分析**：分析基因表达数据集以找到相关途径
- **物种比较**：比较不同生物体的途径数据

**API 基本 URL：** `https://reactome.org/AnalysisService`

### 3.reactome2py Python 包

Python 客户端库封装了 Reactome API 调用，以便更轻松地进行编程访问。

**安装：**
```bash
uv pip install reactome2py
```

**注意：** Reactome2py 软件包（版本 3.0.0，2021 年 1 月发布）可以正常运行，但没有得到积极维护。对于最新功能，请考虑使用直接 REST API 调用。

## 查询Pathway数据

### 使用内容服务 REST API

内容服务使用 REST 协议并以 JSON 或纯文本格式返回数据。

**获取数据库版本：**
<<<代码块_1>>>

**查询特定实体：**
<<<代码块_2>>>

**获取途径中的参与分子：**
<<<代码块_3>>>

### 使用reactome2py包

<<<代码块_4>>>

**详细的API端点和参数**，请参考本技能中的`references/api_reference.md`。

## 执行通路分析

### 过度代表性分析

提交基因/蛋白质标识符列表以查找富集途径。

**使用 REST API：**
<<<代码块_5>>>

**按令牌检索分析：**
<<<代码块_6>>>

### 表达数据分析

使用定量值分析基因表达数据集。

**输入格式（TSV，标头以 # 开头）：**
```
#Gene	Sample1	Sample2	Sample3
TP53	2.5	3.1	2.8
BRCA1	1.2	1.5	1.3
EGFR	4.5	4.2	4.8
```

**提交表达数据：**
```python
import requests

# Read TSV file
with open("expression_data.tsv", "r") as f:
    data = f.read()

response = requests.post(
    "https://reactome.org/AnalysisService/identifiers/",
    headers={"Content-Type": "text/plain"},
    data=data
)

result = response.json()
```

### 物种预测

仅使用 `/projection/` 端点将标识符映射到人类路径：

```python
response = requests.post(
    "https://reactome.org/AnalysisService/identifiers/projection/",
    headers={"Content-Type": "text/plain"},
    data=data
)
```

## 可视化结果

通过使用分析令牌构建 URL，可以在 Reactome Pathway Browser 中可视化分析结果：

```python
token = result["summary"]["token"]
pathway_id = "R-HSA-69278"
url = f"https://reactome.org/PathwayBrowser/#{pathway_id}&DTAB=AN&ANALYSIS={token}"
print(f"View results: {url}")
```

## 使用分析令牌

- 分析令牌的有效期为 **7 天**
- 令牌允许检索先前计算的结果而无需重新提交
- 存储令牌以跨会话访问结果
- 使用 `GET /token/{TOKEN}` 端点检索结果

## 数据格式和标识符

### 支持的标识符类型

Reactome 接受各种标识符格式：
- UniProt 种质（例如 P04637）
- 基因符号（例如TP53）
- Ensembl ID（例如 ENSG00000141510）
- EntrezGene ID（例如 7157）
- 小分子 ChEBI ID

系统自动检测标识符类型。

### 输入格式要求

**对于过度代表性分析：**
- 纯文本标识符列表（每行一个）
- OR TSV 格式的单列

**对于表达分析：**
- TSV 格式，强制标题行以“#”开头
- 第 1 列：标识符
- 第 2+ 列：数值表达式值
- 使用句点 (.) 作为小数点分隔符

### 输出格式
所有 API 响应都会返回包含以下内容的 JSON：
- `pathways`：具有统计指标的丰富路径数组
- `summary`：分析元数据和标记
- `entities`：匹配和未映射的标识符
- 统计值：pValue、FDR（错误发现率）

## 帮助脚本

此技能包括 `scripts/reactome_query.py`，一个用于常见 Reactome 操作的帮助程序脚本：

```bash
# Query pathway information
python scripts/reactome_query.py query R-HSA-69278

# Perform overrepresentation analysis
python scripts/reactome_query.py analyze gene_list.txt

# Get database version
python scripts/reactome_query.py version
```

## 其他资源

- **API 文档**：https://reactome.org/dev
- **用户指南**：https://reactome.org/userguide
- **文档门户**：https://reactome.org/documentation
- **数据下载**：https://reactome.org/download-data
- **reactome2py 文档**：https://reactome.github.io/reactome2py/

有关全面的 API 端点文档，请参阅本技能中的 `references/api_reference.md`。

## 当前数据库统计信息（版本 94，2025 年 9 月）

- 2,825 条人类路径
- 16,002 条反应
- 11,630 种蛋白质
- 2,176 个小分子
- 1,070 种药物
- 41,373 篇参考文献