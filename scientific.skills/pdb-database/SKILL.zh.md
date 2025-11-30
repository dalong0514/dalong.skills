<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pdb数据库
描述：“访问 RCSB PDB 以获取 3D 蛋白质/核酸结构。按文本/序列/结构搜索、下载坐标 (PDB/mmCIF)、检索元数据，以进行结构生物学和药物发现。”
---

# PDB 数据库

## 概述

RCSB PDB 是生物大分子 3D 结构数据的全球存储库。搜索结构、检索坐标和元数据，在 200,000 多个实验确定的结构和计算模型中执行序列和结构相似性搜索。

## 何时使用此技能

该技能应该在以下情况下使用：
- 通过文本、序列或结构相似性搜索蛋白质或核酸 3D 结构
- 下载 PDB、mmCIF 或 BinaryCIF 格式的坐标文件
- 检索结构元数据、实验方法或质量指标
- 跨多个结构执行批量操作
- 将 PDB 数据集成到药物发现、蛋白质工程或结构生物学研究的计算工作流程中

## 核心能力

### 1. 搜索结构

使用各种搜索条件查找 PDB 条目：

**文本搜索：** 按蛋白质名称、关键字或描述搜索
```python
from rcsbapi.search import TextQuery
query = TextQuery("hemoglobin")
results = list(query())
print(f"Found {len(results)} structures")
```

**属性搜索：**查询特定属性（有机体、分辨率、方法等）
<<<代码块_1>>>

**序列相似性：** 查找与给定序列相似的结构
<<<代码块_2>>>

**结构相似性：** 查找具有相似 3D 几何形状的结构
<<<代码块_3>>>

**组合查询：** 使用逻辑运算符构建复杂的搜索
<<<代码块_4>>>

### 2. 检索结构数据

访问有关特定 PDB 条目的详细信息：

**基本入境信息：**
<<<代码块_5>>>

**聚合物实体信息：**
<<<代码块_6>>>

**使用 GraphQL 进行灵活查询：**
```python
from rcsbapi.data import fetch

# Custom GraphQL query
query = """
{
  entry(entry_id: "4HHB") {
    struct {
      title
    }
    exptl {
      method
    }
    rcsb_entry_info {
      resolution_combined
      deposited_atom_count
    }
  }
}
"""
data = fetch(query_type="graphql", query=query)
```

### 3.下载结构文件

检索各种格式的坐标文件：

**下载方法：**
- **PDB 格式**（旧文本格式）：`https://files.rcsb.org/download/{PDB_ID}.pdb`
- **mmCIF 格式**（现代标准）：`https://files.rcsb.org/download/{PDB_ID}.cif`
- **BinaryCIF**（压缩二进制）：使用ModelServer API进行高效访问
- **生物组装**：`https://files.rcsb.org/download/{PDB_ID}.pdb1`（对于组装 1）

**示例下载：**
```python
import requests

pdb_id = "4HHB"

# Download PDB format
pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
response = requests.get(pdb_url)
with open(f"{pdb_id}.pdb", "w") as f:
    f.write(response.text)

# Download mmCIF format
cif_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
response = requests.get(cif_url)
with open(f"{pdb_id}.cif", "w") as f:
    f.write(response.text)
```

### 4. 使用结构数据

检索结构的常见操作：

**解析和分析坐标：**
使用 BioPython 或其他结构生物学库来处理下载的文件：
```python
from Bio.PDB import PDBParser

parser = PDBParser()
structure = parser.get_structure("protein", "4HHB.pdb")

# Iterate through atoms
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(atom.get_coord())
```

**提取元数据：**
```python
from rcsbapi.data import fetch, Schema

# Get experimental details
data = fetch("4HHB", schema=Schema.ENTRY)

resolution = data.get("rcsb_entry_info", {}).get("resolution_combined")
method = data.get("exptl", [{}])[0].get("method")
deposition_date = data.get("rcsb_accession_info", {}).get("deposit_date")

print(f"Resolution: {resolution} Å")
print(f"Method: {method}")
print(f"Deposited: {deposition_date}")
```

### 5. 批量操作

高效处理多种结构：

```python
from rcsbapi.data import fetch, Schema

pdb_ids = ["4HHB", "1MBN", "1GZX"]  # Hemoglobin, myoglobin, etc.

results = {}
for pdb_id in pdb_ids:
    try:
        data = fetch(pdb_id, schema=Schema.ENTRY)
        results[pdb_id] = {
            "title": data["struct"]["title"],
            "resolution": data.get("rcsb_entry_info", {}).get("resolution_combined"),
            "organism": data.get("rcsb_entity_source_organism", [{}])[0].get("scientific_name")
        }
    except Exception as e:
        print(f"Error fetching {pdb_id}: {e}")

# Display results
for pdb_id, info in results.items():
    print(f"\n{pdb_id}: {info['title']}")
    print(f"  Resolution: {info['resolution']} Å")
    print(f"  Organism: {info['organism']}")
```

## Python 包安装

安装官方 RCSB PDB Python API 客户端：

```bash
# Current recommended package
uv pip install rcsb-api

# For legacy code (deprecated, use rcsb-api instead)
uv pip install rcsbsearchapi
```

`rcsb-api` 包通过 `rcsbapi.search` 和 `rcsbapi.data` 模块提供对搜索和数据 API 的统一访问。

## 常见用例

### 药物发现
- 寻找药物靶点的结构
- 分析配体结合位点
- 比较蛋白质-配体复合物
- 识别相似的装订口袋

### 蛋白质工程
- 寻找同源结构进行建模
- 分析序列结构关系
- 比较突变结构
- 研究蛋白质稳定性和动力学

### 结构生物学研究
- 下载结构进行计算分析
- 建立基于结构的对齐方式
- 分析结构特征（二级结构、域）
- 比较实验方法和质量指标

### 教育和可视化
- 检索教学结构
- 生成分子可视化
- 探索结构与功能的关系
- 研究进化保护

## 关键概念

**PDB ID：** 每个结构条目的唯一 4 字符标识符（例如“4HHB”）。 AlphaFold 和 ModelArchive 条目以“AF_”或“MA_”前缀开头。

**mmCIF/PDBx：** 使用键值结构的现代文件格式，替换大型结构的旧 PDB 格式。

**生物组装：** 大分子的功能形式，可能包含来自不对称单元的多个链副本。

**分辨率：** 晶体结构细节的测量（较低的值=较高的细节）。典型范围：1.5-3.5 Å，用于高质量结构。

**实体：** 结构中独特的分子成分（蛋白质链、DNA、配体等）。

## 资源

此技能包括 `references/` 目录中的参考文档：

### 参考资料/api_reference.md
全面的 API 文档涵盖：
- 详细的API端点规范
- 高级查询模式和示例
- 数据模式参考
- 速率限制和最佳实践
- 解决常见问题

当您需要有关 API 功能、复杂查询构造或详细数据架构信息的深入信息时，请使用此参考。

## 其他资源

- **RCSB PDB 网站：** https://www.rcsb.org
- **PDB-101 教育门户：** https://pdb101.rcsb.org
- **API 文档：** https://www.rcsb.org/docs/programmatic-access/web-apis-overview
- **Python 包文档：** https://rcsbapi.readthedocs.io/
- **数据 API 文档：** https://data.rcsb.org/
- **GitHub 存储库：** https://github.com/rcsb/py-rcsb-api