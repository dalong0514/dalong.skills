<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 锌数据库
描述：“访问 ZINC（超过 2.3 亿可购买化合物）。通过 ZINC ID/SMILES 搜索、相似性搜索、用于对接、类似物发现、虚拟筛选和药物发现的 3D 就绪结构。”
---

# ZINC 数据库

## 概述

ZINC 是一个可免费访问的存储库，包含 UCSF 维护的 2.3 亿多个可购买化合物。通过 ZINC ID 或 SMILES 进行搜索，执行相似性搜索，下载 3D 就绪结构进行对接，发现用于虚拟筛选和药物发现的类似物。

## 何时使用此技能

该技能应该在以下情况下使用：

- **虚拟筛选**：寻找用于分子对接研究的化合物
- **先导化合物发现**：识别用于药物开发的商用化合物
- **结构搜索**：通过 SMILES 执行相似性或类似搜索
- **化合物检索**：通过 ZINC ID 或供应商代码查找分子
- **化学空间探索**：探索可购买的化学多样性
- **对接研究**：访问 3D 就绪分子结构
- **模拟搜索**：根据结构相似性查找相似化合物
- **供应商查询**：识别来自特定化学品供应商的化合物
- **随机抽样**：获取随机化合物集进行筛选

## 数据库版本

ZINC 已经发展了多个版本：

- **ZINC22**（当前）：最大版本，拥有超过 230 百万种可购买化合物和数十亿规模的按需生产化合物
- **ZINC20**：仍然保留，专注于类铅和类药物化合物
- **ZINC15**：前身版本，遗留但仍记录在案

该技能主要关注 ZINC22，这是最新、最全面的版本。

## 访问方法

### 网页界面

主要接入点：https://zinc.docking.org/
交互式搜索：https://cartblanche22.docking.org/

### API 访问

所有 ZINC22 搜索都可以通过 CartBlanche22 API 以编程方式执行：

**基本网址**：`https://cartblanche22.docking.org/`

所有 API 端点都以文本或 JSON 格式返回带有可自定义字段的数据。

## 核心能力

### 1. 按 ZINC ID 搜索

使用 ZINC 标识符检索特定化合物。

**网络界面**：https://cartblanche22.docking.org/search/zincid

**API端点**：
```bash
curl "https://cartblanche22.docking.org/[email protected]_fields=smiles,zinc_id"
```

**多个ID**：
<<<代码块_1>>>

**响应字段**：`zinc_id`、`smiles`、`sub_id`、`supplier_code`、`catalogs`、`tranche`（包括 H-count、LogP、MW、相位）

### 2. 按微笑搜索

使用 SMILES 表示法按化学结构查找化合物，并具有用于模拟搜索的可选距离参数。

**网络界面**：https://cartblanche22.docking.org/search/smiles

**API端点**：
<<<代码块_2>>>

**参数**：
- `smiles`：查询 SMILES 字符串（如有必要，进行 URL 编码）
- `dist`：谷本距离阈值（默认值：0 表示完全匹配）
- `adist`：用于更广泛搜索的替代距离参数（默认值：0）
- `output_fields`：所需输出字段的逗号分隔列表

**示例 - 完全匹配**：
<<<代码块_3>>>

**示例 - 相似性搜索**：
<<<代码块_4>>>

### 3. 按供应商代码搜索

查询来自特定化学品供应商的化合物或从特定目录中检索所有分子。

**网页界面**：https://cartblanche22.docking.org/search/catitems

**API端点**：
<<<代码块_5>>>

**用例**：
- 验证特定供应商的化合物可用性
- 从目录中检索所有化合物
- 交叉引用供应商代码与 ZINC ID

### 4. 随机复合采样

生成随机化合物集用于筛选或基准测试目的。

**网页界面**：https://cartblanche22.docking.org/search/random

**API端点**：
<<<代码块_6>>>

**参数**：
- `count`：要检索的随机化合物数（默认值：100）
- `subset`：按子集过滤（例如，“类铅”、“类药物”、“片段”）
- `output_fields`：自定义返回的数据字段

**示例 - 随机铅状分子**：
```bash
curl "https://cartblanche22.docking.org/substance/random.txt:count=1000&subset=lead-like&output_fields=zinc_id,smiles,tranche"
```

## 常见工作流程

### 工作流程 1：准备对接库

1. **根据目标特性或所需的化学空间定义搜索标准**

2. **使用适当的搜索方法查询 ZINC22**：
   ```bash
   # Example: Get drug-like compounds with specific LogP and MW
   curl "https://cartblanche22.docking.org/substance/random.txt:count=10000&subset=drug-like&output_fields=zinc_id,smiles,tranche" > docking_library.txt
   ```

3. **解析结果**以提取 ZINC ID 和 SMILES：
   ```python
   import pandas as pd

   # Load results
   df = pd.read_csv('docking_library.txt', sep='\t')

   # Filter by properties in tranche data
   # Tranche format: H##P###M###-phase
   # H = H-bond donors, P = LogP*10, M = MW
   ```

4. **使用 ZINC ID 下载 3D 结构**以进行对接或从文件存储库下载

### 工作流程 2：寻找热门化合物的类似物

1. **获取命中化合物的 SMILES**：
   ```python
   hit_smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Example: Ibuprofen
   ```

2. **使用距离阈值执行相似性搜索**：
   ```bash
   curl "https://cartblanche22.docking.org/smiles.txt:smiles=CC(C)Cc1ccc(cc1)C(C)C(=O)O&dist=5&output_fields=zinc_id,smiles,catalogs" > analogs.txt
   ```
3. **分析结果**以识别可购买的类似物：
   ```python
   import pandas as pd

   analogs = pd.read_csv('analogs.txt', sep='\t')
   print(f"Found {len(analogs)} analogs")
   print(analogs[['zinc_id', 'smiles', 'catalogs']].head(10))
   ```

4. **检索 3D 结构**以获得最有前途的类似物

### 工作流程 3：批量化合物检索

1. **从文献、数据库或之前的屏幕中编译 ZINC ID 列表**：
   ```python
   zinc_ids = [
       "ZINC000000000001",
       "ZINC000000000002",
       "ZINC000000000003"
   ]
   zinc_ids_str = ",".join(zinc_ids)
   ```

2. **查询ZINC22 API**：
   ```bash
   curl "https://cartblanche22.docking.org/substances.txt:zinc_id=ZINC000000000001,ZINC000000000002&output_fields=zinc_id,smiles,supplier_code,catalogs"
   ```

3. **处理结果**用于下游分析或采购

### 工作流程 4：化学空间采样

1. **根据筛选目标选择子集参数**：
   - 片段：MW < 250，有利于基于片段的药物发现
   - 类铅：MW 250-350，LogP ≤ 3.5
   - 类药物：MW 350-500，遵循利平斯基五法则

2. **生成随机样本**：
   ```bash
   curl "https://cartblanche22.docking.org/substance/random.txt:count=5000&subset=lead-like&output_fields=zinc_id,smiles,tranche" > chemical_space_sample.txt
   ```

3. **分析化学多样性**并为虚拟筛选做好准备

## 输出字段

使用 `output_fields` 参数自定义 API 响应：

**可用字段**：
- `zinc_id`：ZINC 标识符
- `smiles`：SMILES 字符串表示形式
- `sub_id`：内部物质 ID
- `supplier_code`：供应商目录号
- `catalogs`：提供化合物的供应商列表
- `tranche`：编码分子特性（H 计数、LogP、MW、反应相）

**示例**：
```bash
curl "https://cartblanche22.docking.org/substances.txt:zinc_id=ZINC000000000001&output_fields=zinc_id,smiles,catalogs,tranche"
```

## 分级制度

ZINC 根据分子特性将化合物组织成“部分”：

**格式**：`H##P###M###-phase`

- **H##**：氢键供体数量 (00-99)
- **P###**：LogP × 10（例如，P035 = LogP 3.5）
- **M###**：以道尔顿为单位的分子量（例如，M400 = 400 Da）
- **阶段**：反应性分类

**部分示例**：`H05P035M400-0`
- 5 个氢键捐助者
- LogP = 3.5
- 兆瓦 = 400 Da
- 反应阶段 0

使用批次数据按药物相似性标准过滤化合物。

## 下载 3D 结构

对于分子对接，可通过文件存储库获取 3D 结构：

**文件存储库**：https://files.docking.org/zinc22/

结构按部分组织，并以多种格式提供：
- MOL2：具有 3D 坐标的多分子格式
- SDF：结构数据文件格式
- DB2.GZ：DOCK 的压缩数据库格式

请参阅 https://wiki.docking.org 处的 ZINC 文档，了解下载协议和批量访问方法。

## Python 集成

### 在 Python 中使用curl

```python
import subprocess
import json

def query_zinc_by_id(zinc_id, output_fields="zinc_id,smiles,catalogs"):
    """Query ZINC22 by ZINC ID."""
    url = f"https://cartblanche22.docking.org/[email protected]_id={zinc_id}&output_fields={output_fields}"
    result = subprocess.run(['curl', url], capture_output=True, text=True)
    return result.stdout

def search_by_smiles(smiles, dist=0, adist=0, output_fields="zinc_id,smiles"):
    """Search ZINC22 by SMILES with optional distance parameters."""
    url = f"https://cartblanche22.docking.org/smiles.txt:smiles={smiles}&dist={dist}&adist={adist}&output_fields={output_fields}"
    result = subprocess.run(['curl', url], capture_output=True, text=True)
    return result.stdout

def get_random_compounds(count=100, subset=None, output_fields="zinc_id,smiles,tranche"):
    """Get random compounds from ZINC22."""
    url = f"https://cartblanche22.docking.org/substance/random.txt:count={count}&output_fields={output_fields}"
    if subset:
        url += f"&subset={subset}"
    result = subprocess.run(['curl', url], capture_output=True, text=True)
    return result.stdout
```

### 解析结果

```python
import pandas as pd
from io import StringIO

# Query ZINC and parse as DataFrame
result = query_zinc_by_id("ZINC000000000001")
df = pd.read_csv(StringIO(result), sep='\t')

# Extract tranche properties
def parse_tranche(tranche_str):
    """Parse ZINC tranche code to extract properties."""
    # Format: H##P###M###-phase
    import re
    match = re.match(r'H(\d+)P(\d+)M(\d+)-(\d+)', tranche_str)
    if match:
        return {
            'h_donors': int(match.group(1)),
            'logP': int(match.group(2)) / 10.0,
            'mw': int(match.group(3)),
            'phase': int(match.group(4))
        }
    return None

df['tranche_props'] = df['tranche'].apply(parse_tranche)
```

## 最佳实践

### 查询优化

- **开始具体**：从精确搜索开始，然后扩展到相似性搜索
- **使用适当的距离参数**：较小的距离值 (1-3) 用于接近的类似物，较大的距离值 (5-10) 用于不同的类似物
- **限制输出字段**：仅请求必要的字段以减少数据传输
- **批量查询**：如果可能，在单个 API 调用中组合多个 ZINC ID

### 性能考虑因素

- **速率限制**：尊重服务器资源；避免快速连续的请求
- **缓存**：在本地存储经常访问的化合物
- **并行下载**：下载 3D 结构时，对文件存储库使用并行 wget 或 aria2c
- **子集过滤**：使用铅样、药物样或片段子集来减少搜索空间

### 数据质量

- **验证可用性**：供应商目录发生变化；在大订单之前确认化合物的可用性
- **检查立体化学**：SMILES 可能无法完全指定立体化学；验证 3D 结构
- **验证结构**：使用化学信息学工具（RDKit、OpenBabel）验证结构有效性
- **交叉参考**：如果可能，与其他数据库（PubChem、ChEMBL）交叉检查

## 资源

### 参考文献/api_reference.md

综合文档包括：

- 完整的API端点参考
- URL语法和参数规范
- 高级查询模式和示例
- 文件存储库组织和访问
- 批量下载方法
- 错误处理和故障排除
- 与分子对接软件集成

请参阅本文档以获取详细的技术信息和高级使用模式。

## 重要免责声明

### 数据可靠性

ZINC 明确声明：**“我们不保证任何用途的任何分子的质量，也不对因使用该数据库而产生的错误承担任何责任。”**

- 化合物的供应情况可能会发生变化，恕不另行通知
- 结构表示可能包含错误
- 供应商信息应独立验证
- 在实验工作之前进行适当的验证

### 适当使用

- ZINC 用于药物发现的学术和研究目的
- 验证商业用途的许可条款
- 使用专利化合物时尊重知识产权
- 遵循您所在机构的化合物采购指南

## 其他资源

- **ZINC 网站**：https://zinc.docking.org/
- **CartBlanche22 界面**：https://cartblanche22.docking.org/
- **ZINC 维基**：https://wiki.docking.org/
- **文件存储库**：https://files.docking.org/zinc22/
- **GitHub**：https://github.com/docking-org/
- **主要出版物**：Irwin 等人，J. Chem。信息。 2020 型 (ZINC15)
- **ZINC22 出版物**：Irwin 等人，J. Chem。信息。型号2023

## 引文

在出版物中使用 ZINC 时，请引用适当的版本：

**锌22**：
欧文，J.J.，等人。 “ZINC22——用于配体发现的免费的数十亿规模的有形化合物数据库。” *化学信息与建模杂志* 2023。

**锌15**：
欧文，J.J.，等人。 “ZINC15——适合所有人的配体发现。” *化学信息与建模杂志* 2020, 60, 6065–6073。