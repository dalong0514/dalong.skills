<!-- 此文件由机器翻译自 api_reference.md -->

# ZINC 数据库 API 参考

## 概述

以编程方式访问 ZINC 数据库的完整技术参考，涵盖 ZINC22、ZINC20 和旧版本的 API 端点、查询语法、参数、响应格式以及高级使用模式。

## 基本 URL

### ZINC22（当前）
- **CartBlanche22 API**：`https://cartblanche22.docking.org/`
- **文件存储库**：`https://files.docking.org/zinc22/`
- **主网站**：`https://zinc.docking.org/`

### ZINC20（维持）
- **API**：`https://zinc20.docking.org/`
- **文件存储库**：`https://files.docking.org/zinc20/`

### 文档
- **维基**：`https://wiki.docking.org/`
- **GitHub**：`https://github.com/docking-org/`

## API 端点

### 1. 通过 ZINC ID 检索物质

使用 ZINC 标识符检索化合物信息。

**端点**：`/substances.txt`

**参数**：
- `zinc_id`（必需）：单个 ZINC ID 或逗号分隔列表
- `output_fields`（可选）：逗号分隔的字段名称（默认值：所有字段）

**网址格式**：
```
https://cartblanche22.docking.org/substances.txt:zinc_id={ZINC_ID}&output_fields={FIELDS}
```

**示例**：

单一化合物：
<<<代码块_1>>>

多种化合物：
<<<代码块_2>>>

从文件中批量检索：
<<<代码块_3>>>

**响应格式** (TSV)：
<<<代码块_4>>>

### 2. SMILES 结构搜索

通过具有可选相似性阈值的化学结构搜索化合物。

**端点**：`/smiles.txt`

**参数**：
- `smiles`（必需）：查询 SMILES 字符串（必要时进行 URL 编码）
- `dist`（可选）：Tanimoto 距离阈值（0-10，默认值：0 = 精确）
- `adist`（可选）：替代距离度量（0-10，默认值：0）
- `output_fields`（可选）：逗号分隔的字段名称

**网址格式**：
<<<代码块_5>>>

**示例**：

结构精确匹配：
<<<代码块_6>>>

相似性搜索（谷本距离 = 3）：
```bash
curl "https://cartblanche22.docking.org/smiles.txt:smiles=CC(C)Cc1ccc(cc1)C(C)C(=O)O&dist=3&output_fields=zinc_id,smiles,catalogs"
```

广泛的相似性搜索：
```bash
curl "https://cartblanche22.docking.org/smiles.txt:smiles=c1ccccc1&dist=5&adist=5&output_fields=zinc_id,smiles,tranche"
```

URL 编码的 SMILES（对于特殊字符）：
```bash
# Original: CC(=O)Oc1ccccc1C(=O)O
# Encoded: CC%28%3DO%29Oc1ccccc1C%28%3DO%29O
curl "https://cartblanche22.docking.org/smiles.txt:smiles=CC%28%3DO%29Oc1ccccc1C%28%3DO%29O&dist=2"
```

**距离参数解释**：
- `dist=0`：完全匹配
- `dist=1-3`：相近的类似物（高度相似性）
- `dist=4-6`：中等类似物
- `dist=7-10`：多样化的化学空间

### 3. 供应商代码搜索

按供应商目录号查询化合物。

**端点**：`/catitems.txt`

**参数**：
- `catitem_id`（必需）：供应商目录代码
- `output_fields`（可选）：逗号分隔的字段名称

**网址格式**：
```
https://cartblanche22.docking.org/catitems.txt:catitem_id={SUPPLIER_CODE}&output_fields={FIELDS}
```

**示例**：
```bash
curl "https://cartblanche22.docking.org/catitems.txt:catitem_id=SUPPLIER-12345&output_fields=zinc_id,smiles,supplier_code,catalogs"
```

### 4. 随机复合采样

生成随机化合物集，并可选择按化学性质进行过滤。

**端点**：`/substance/random.txt`

**参数**：
- `count`（可选）：要检索的化合物数量（默认值：100，最大值：取决于服务器）
- `subset`（可选）：按预定义子集过滤（例如，“类先导”、“类药物”、“片段”）
- `output_fields`（可选）：逗号分隔的字段名称

**网址格式**：
```
https://cartblanche22.docking.org/substance/random.txt:count={COUNT}&subset={SUBSET}&output_fields={FIELDS}
```

**示例**：

随机 100 个化合物（默认）：
```bash
curl "https://cartblanche22.docking.org/substance/random.txt"
```

随机铅样分子：
```bash
curl "https://cartblanche22.docking.org/substance/random.txt:count=1000&subset=lead-like&output_fields=zinc_id,smiles,tranche"
```

随机药物样分子：
```bash
curl "https://cartblanche22.docking.org/substance/random.txt:count=5000&subset=drug-like&output_fields=zinc_id,smiles"
```

随机片段：
```bash
curl "https://cartblanche22.docking.org/substance/random.txt:count=500&subset=fragment&output_fields=zinc_id,smiles,tranche"
```

**子集定义**：
- `fragment`：MW < 250，适合基于片段的药物发现
- `lead-like`：MW 250-350，LogP ≤ 3.5，可旋转键 ≤ 7
- `drug-like`：MW 350-500，遵循 Lipinski 的五法则
- `lugs`：大型且异常优秀的子集（精心策划）

## 输出字段

### 可用字段

使用 `output_fields` 参数自定义 API 响应：

|领域 |描述 |示例|
|--------|-------------|---------|
| `zinc_id` |锌标识符| ZINC000000000001 |
| `smiles` |规范 SMILES 字符串 | CC(C)O |
| `sub_id` |内部物质 ID | 123456 |
| `supplier_code` |供应商目录号 | AB-1234567 |
| `catalogs` |供应商名录| [分子、mcule、mcule-ultimate] |
| `tranche` |编码分子特性 | H02P025M300-0 |
| `mwt` |分子量| 325.45 | 325.45
| `logp` | LogP（分配系数）| 2.5 | 2.5
| `hba` | H键受体| 4 |
| `hbd` | H 型债券捐助者 | 2 |
| `rotatable_bonds` |可旋转债券数量 | 5 |

**注意**：并非所有字段都适用于所有端点。字段可用性取决于数据库版本和端点。

### 默认字段
如果未指定 `output_fields`，端点将返回 TSV 格式的所有可用字段。

### 自定义字段选择

仅请求特定字段：
```bash
curl "https://cartblanche22.docking.org/[email protected]_fields=zinc_id,smiles"
```

请求多个字段：
```bash
curl "https://cartblanche22.docking.org/[email protected]_fields=zinc_id,smiles,tranche,catalogs"
```

## 分级制度

ZINC 根据分子特性将化合物组织成级，以实现高效过滤和组织。

### 批次代码格式

**模式**：`H##P###M###-phase`

|组件|描述 |范围 |
|------------|-------------|--------|
| H## |氢键供体| 00-99|
| P### | LogP × 10 | 000-999（例如，P035 = LogP 3.5）|
| M### |分子量| 000-999 大 |
|相|反应性分类| 0-9 |

### 示例

|批次代码 |解读|
|--------------|----------------|
| `H00P010M250-0` | 0 H-供体，LogP=1.0，MW=250 Da，阶段 0 |
| `H05P035M400-0` | 5 个 H 供体，LogP=3.5，MW=400 Da，阶段 0 |
| `H02P-005M180-0` | 2 个 H 供体，LogP=-0.5，MW=180 Da，阶段 0 |

### 反应阶段

|相|描述 |
|--------|-------------|
| 0 |无反应（首选筛查）|
| 1-9 | 1-9增加反应性（PAINS、反应基团）|

### 在 Python 中解析 Tranches

```python
import re

def parse_tranche(tranche_str):
    """
    Parse ZINC tranche code.

    Args:
        tranche_str: Tranche code (e.g., "H05P035M400-0")

    Returns:
        dict with h_donors, logp, mw, phase
    """
    pattern = r'H(\d+)P(-?\d+)M(\d+)-(\d+)'
    match = re.match(pattern, tranche_str)

    if not match:
        return None

    return {
        'h_donors': int(match.group(1)),
        'logp': int(match.group(2)) / 10.0,
        'mw': int(match.group(3)),
        'phase': int(match.group(4))
    }

# Example usage
tranche = "H05P035M400-0"
props = parse_tranche(tranche)
print(props)  # {'h_donors': 5, 'logp': 3.5, 'mw': 400, 'phase': 0}
```

### 按批次过滤

从文件存储库下载特定部分：
```bash
# Download all compounds in a specific tranche
wget https://files.docking.org/zinc22/H05/H05P035M400-0.db2.gz
```

## 文件存储库访问

### 目录结构

ZINC22 3D 结构按氢键供体分层组织：

```
https://files.docking.org/zinc22/
├── H00/
│   ├── H00P010M200-0.db2.gz
│   ├── H00P020M250-0.db2.gz
│   └── ...
├── H01/
├── H02/
└── ...
```

### 文件格式

|扩展|格式|描述 |
|------------|--------|-------------|
| `.db2.gz` | DOCK数据库| DOCK 的压缩多一致性数据库 |
| `.mol2.gz` |摩尔2 |具有 3D 坐标的多分子格式 |
| `.sdf.gz` |自卫队|结构-数据文件格式|
| `.smi` |微笑|带有 ZINC ID 的纯文本 SMILES |

### 下载 3D 结构

**单一批次**：
```bash
wget https://files.docking.org/zinc22/H05/H05P035M400-0.db2.gz
```

**多部分**（与 aria2c 并行下载）：
```bash
# Create URL list
cat > tranche_urls.txt <<EOF
https://files.docking.org/zinc22/H05/H05P035M400-0.db2.gz
https://files.docking.org/zinc22/H05/H05P035M400-0.db2.gz
https://files.docking.org/zinc22/H05/H05P040M400-0.db2.gz
EOF

# Download in parallel
aria2c -i tranche_urls.txt -x 8 -j 4
```

**递归下载**（谨慎使用-大数据）：
```bash
wget -r -np -nH --cut-dirs=1 -A "*.db2.gz" \
  https://files.docking.org/zinc22/H05/
```

### 提取结构

```bash
# Decompress
gunzip H05P035M400-0.db2.gz

# Convert to other formats using OpenBabel
obabel H05P035M400-0.db2 -O output.sdf
obabel H05P035M400-0.db2 -O output.mol2
```

## 高级查询模式

### 组合多个搜索条件

**复杂查询的 Python 包装器**：

```python
import subprocess
import pandas as pd
from io import StringIO

def advanced_zinc_search(smiles=None, zinc_ids=None, dist=0,
                         subset=None, count=None, output_fields=None):
    """
    Flexible ZINC search with multiple criteria.

    Args:
        smiles: SMILES string for structure search
        zinc_ids: List of ZINC IDs for batch retrieval
        dist: Distance parameter for similarity (0-10)
        subset: Subset filter (lead-like, drug-like, fragment)
        count: Number of random compounds
        output_fields: List of fields to return

    Returns:
        pandas DataFrame with results
    """
    if output_fields is None:
        output_fields = ['zinc_id', 'smiles', 'tranche', 'catalogs']

    fields_str = ','.join(output_fields)

    # Structure search
    if smiles:
        url = f"https://cartblanche22.docking.org/smiles.txt:smiles={smiles}&dist={dist}&output_fields={fields_str}"

    # Batch retrieval
    elif zinc_ids:
        zinc_ids_str = ','.join(zinc_ids)
        url = f"https://cartblanche22.docking.org/substances.txt:zinc_id={zinc_ids_str}&output_fields={fields_str}"

    # Random sampling
    elif count:
        url = f"https://cartblanche22.docking.org/substance/random.txt:count={count}&output_fields={fields_str}"
        if subset:
            url += f"&subset={subset}"

    else:
        raise ValueError("Must specify smiles, zinc_ids, or count")

    # Execute query
    result = subprocess.run(['curl', '-s', url],
                          capture_output=True, text=True)

    # Parse to DataFrame
    df = pd.read_csv(StringIO(result.stdout), sep='\t')

    return df
```

**使用示例**：

```python
# Find similar compounds
df = advanced_zinc_search(
    smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    dist=3,
    output_fields=['zinc_id', 'smiles', 'catalogs']
)

# Batch retrieval
zinc_ids = ["ZINC000000000001", "ZINC000000000002"]
df = advanced_zinc_search(zinc_ids=zinc_ids)

# Random drug-like set
df = advanced_zinc_search(
    count=1000,
    subset='drug-like',
    output_fields=['zinc_id', 'smiles', 'tranche']
)
```

### 基于属性的过滤

使用批次数据按分子特性过滤化合物：

```python
def filter_by_properties(df, mw_range=None, logp_range=None,
                        max_hbd=None, phase=0):
    """
    Filter DataFrame by molecular properties.

    Args:
        df: DataFrame with 'tranche' column
        mw_range: Tuple (min_mw, max_mw)
        logp_range: Tuple (min_logp, max_logp)
        max_hbd: Maximum H-bond donors
        phase: Reactivity phase (0 = unreactive)

    Returns:
        Filtered DataFrame
    """
    # Parse tranches
    df['tranche_props'] = df['tranche'].apply(parse_tranche)
    df['mw'] = df['tranche_props'].apply(lambda x: x['mw'] if x else None)
    df['logp'] = df['tranche_props'].apply(lambda x: x['logp'] if x else None)
    df['hbd'] = df['tranche_props'].apply(lambda x: x['h_donors'] if x else None)
    df['phase'] = df['tranche_props'].apply(lambda x: x['phase'] if x else None)

    # Apply filters
    mask = pd.Series([True] * len(df))

    if mw_range:
        mask &= (df['mw'] >= mw_range[0]) & (df['mw'] <= mw_range[1])

    if logp_range:
        mask &= (df['logp'] >= logp_range[0]) & (df['logp'] <= logp_range[1])

    if max_hbd is not None:
        mask &= df['hbd'] <= max_hbd

    if phase is not None:
        mask &= df['phase'] == phase

    return df[mask]

# Example: Get drug-like compounds with specific properties
df = advanced_zinc_search(count=10000, subset='drug-like')
filtered = filter_by_properties(
    df,
    mw_range=(300, 450),
    logp_range=(1.0, 4.0),
    max_hbd=3,
    phase=0
)
```

## 速率限制和最佳实践

### 速率限制

ZINC 没有发布明确的速率限制，但用户应该：

- **避免快速请求**：将查询间隔至少 1 秒
- **使用批量操作**：在单个请求中查询多个 ZINC ID
- **缓存结果**：在本地存储经常访问的数据
- **非高峰使用**：在非高峰时段（UTC 夜间/周末）执行大量下载

### 礼仪

```python
import time

def polite_zinc_query(query_func, *args, delay=1.0, **kwargs):
    """Wrapper to add delay between queries."""
    result = query_func(*args, **kwargs)
    time.sleep(delay)
    return result
```

### 错误处理

```python
def robust_zinc_query(url, max_retries=3, timeout=30):
    """
    Query ZINC with retry logic.

    Args:
        url: Full ZINC API URL
        max_retries: Maximum retry attempts
        timeout: Request timeout in seconds

    Returns:
        Query results or None on failure
    """
    import subprocess
    import time

    for attempt in range(max_retries):
        try:
            result = subprocess.run(
                ['curl', '-s', '--max-time', str(timeout), url],
                capture_output=True,
                text=True,
                check=True
            )

            # Check for empty or error responses
            if not result.stdout or 'error' in result.stdout.lower():
                raise ValueError("Invalid response")

            return result.stdout

        except (subprocess.CalledProcessError, ValueError) as e:
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt  # Exponential backoff
                print(f"Retry {attempt + 1}/{max_retries} after {wait_time}s...")
                time.sleep(wait_time)
            else:
                print(f"Failed after {max_retries} attempts")
                return None
```

## 与分子对接整合

### 准备 DOCK6 库

```bash
# 1. Download tranche files
wget https://files.docking.org/zinc22/H05/H05P035M400-0.db2.gz

# 2. Decompress
gunzip H05P035M400-0.db2.gz

# 3. Use directly with DOCK6
dock6 -i dock.in -o dock.out -l H05P035M400-0.db2
```

### AutoDock Vina 集成

```bash
# 1. Download MOL2 format
wget https://files.docking.org/zinc22/H05/H05P035M400-0.mol2.gz
gunzip H05P035M400-0.mol2.gz

# 2. Convert to PDBQT using prepare_ligand script
prepare_ligand4.py -l H05P035M400-0.mol2 -o ligands.pdbqt -A hydrogens

# 3. Run Vina
vina --receptor protein.pdbqt --ligand ligands.pdbqt \
     --center_x 25.0 --center_y 25.0 --center_z 25.0 \
     --size_x 20.0 --size_y 20.0 --size_z 20.0
```

### RDKit 集成

```python
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import pandas as pd

def process_zinc_results(zinc_df):
    """
    Process ZINC results with RDKit.

    Args:
        zinc_df: DataFrame with SMILES column

    Returns:
        DataFrame with calculated properties
    """
    # Convert SMILES to molecules
    zinc_df['mol'] = zinc_df['smiles'].apply(Chem.MolFromSmiles)

    # Calculate properties
    zinc_df['mw'] = zinc_df['mol'].apply(Descriptors.MolWt)
    zinc_df['logp'] = zinc_df['mol'].apply(Descriptors.MolLogP)
    zinc_df['hbd'] = zinc_df['mol'].apply(Descriptors.NumHDonors)
    zinc_df['hba'] = zinc_df['mol'].apply(Descriptors.NumHAcceptors)
    zinc_df['tpsa'] = zinc_df['mol'].apply(Descriptors.TPSA)
    zinc_df['rotatable'] = zinc_df['mol'].apply(Descriptors.NumRotatableBonds)

    # Generate 3D conformers
    for mol in zinc_df['mol']:
        if mol:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

    return zinc_df

# Save to SDF for docking
def save_to_sdf(zinc_df, output_file):
    """Save molecules to SDF file."""
    writer = Chem.SDWriter(output_file)
    for idx, row in zinc_df.iterrows():
        if row['mol']:
            row['mol'].SetProp('ZINC_ID', row['zinc_id'])
            writer.write(row['mol'])
    writer.close()
```

## 故障排除

### 常见问题

**问题**：空或无结果
- **解决方案**：检查 SMILES 语法，验证 ZINC ID 是否存在，尝试更广泛的相似性搜索

**问题**：超时错误
- **解决方案**：减少结果数量，使用批量查询，在非高峰时段尝试

**问题**：SMILES 编码无效
- **解决方案**：对特殊字符进行 URL 编码（在 Python 中使用 `urllib.parse.quote()`）

**问题**：找不到批次文件
- **解决方案**：验证批次代码格式，检查文件存储库结构

### 调试模式

```python
def debug_zinc_query(url):
    """Print query details for debugging."""
    print(f"Query URL: {url}")

    result = subprocess.run(['curl', '-v', url],
                          capture_output=True, text=True)

    print(f"Status: {result.returncode}")
    print(f"Stderr: {result.stderr}")
    print(f"Stdout length: {len(result.stdout)}")
    print(f"First 500 chars:\n{result.stdout[:500]}")

    return result.stdout
```

## 版本差异

### ZINC22 vs ZINC20 vs ZINC15

|特色 |锌22 |锌20 |锌15 |
|--------|--------|--------|--------|
|化合物| 230M+ 可购买 |专注于潜在客户|总计约 750M |
|应用程序接口 |购物车Blanche22 |类似|类似休息 |
|部分 |是的 |是的 |是的 |
| 3D 结构 |是的 |是的 |是的 |
|状态 |当前，不断增长|维护|遗产|

### API 兼容性

大多数查询模式可以跨版本工作，但 URL 有所不同：
- ZINC22：`cartblanche22.docking.org`
- ZINC20：`zinc20.docking.org`
- ZINC15：`zinc15.docking.org`

## 其他资源

- **ZINC 维基**：https://wiki.docking.org/
- **ZINC22 文档**：https://wiki.docking.org/index.php/Category:ZINC22
- **ZINC API 指南**：https://wiki.docking.org/index.php/ZINC_api
- **文件访问指南**：https://wiki.docking.org/index.php/ZINC22:Getting_started
- **出版物**：
  - ZINC22：J. Chem。信息。模型。 2023年
  - ZINC15：J. Chem。信息。模型。 2020, 60, 6065-6073
- **支持**：通过 ZINC 网站或 GitHub 问题联系