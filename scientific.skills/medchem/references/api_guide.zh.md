<!-- 此文件由机器翻译自 api_guide.md -->

# Medchem API 参考

所有 medchem 模块和功能的综合参考。

## 模块：medchem.rules

### 类：规则过滤器

根据多种药物化学规则过滤分子。

**构造函数：**
```python
RuleFilters(rule_list: List[str])
```

**参数：**
- `rule_list`：要应用的规则名称列表。请参阅下面的可用规则。

**方法：**

<<<代码块_1>>>
- `mols`：RDKit分子对象列表
- `n_jobs`：并行作业数量（-1 使用所有核心）
- `progress`：显示进度条
- **返回**：包含每个规则结果的字典

**示例：**
<<<代码块_2>>>

### 模块：medchem.rules.basic_rules

可应用于单个分子的单独规则函数。

#### 五的规则（）

<<<代码块_3>>>

口服生物利用度的利平斯基五法则。

**标准：**
- 分子量 ≤ 500 Da
- LogP ≤ 5
-氢键供体≤5
- H键受体≤10

**参数：**
- `mol`：SMILES 字符串或 RDKit 分子对象

**返回：** 如果分子通过所有标准，则返回 True

#### 三个规则（）

<<<代码块_4>>>

片段筛选文库的三法则。

**标准：**
- 分子量 ≤ 300 Da
- LogP ≤ 3
-氢键供体≤3
- 氢键受体 ≤ 3
- 可旋转键 ≤ 3
- 极性表面积 ≤ 60 Ų

####rule_of_oprea()

<<<代码块_5>>>

Oprea 的类似线索标准，用于点击线索优化。

**标准：**
- 分子量：200-350 Da
- LogP：-2至4
- 可旋转键 ≤ 7
- 戒指 ≤ 4

####rule_of_cns()

<<<代码块_6>>>

CNS 药物相似规则。

**标准：**
- 分子量 ≤ 450 Da
- LogP：-1至5
- 氢键供体 ≤ 2
- TPSA ≤ 90 Ų

####rule_of_leadlike_soft()

```python
rule_of_leadlike_soft(mol: Union[str, Chem.Mol]) -> bool
```

类似软铅的标准（更宽松）。

**标准：**
- 分子量：250-450 Da
- LogP：-3至4
- 可旋转键 ≤ 10

####rule_of_leadlike_strict()

```python
rule_of_leadlike_strict(mol: Union[str, Chem.Mol]) -> bool
```

严格的类似铅标准（更具限制性）。

**标准：**
- 分子量：200-350 Da
- LogP：-2至3.5
- 可旋转键 ≤ 7
- 戒指：1-3

####rule_of_veber()

```python
rule_of_veber(mol: Union[str, Chem.Mol]) -> bool
```

Veber 的口服生物利用度规则。

**标准：**
- 可旋转键 ≤ 10
- TPSA ≤ 140 Ų

####rule_of_reos()

```python
rule_of_reos(mol: Union[str, Chem.Mol]) -> bool
```

快速消除泔水 (REOS) 过滤器。

**标准：**
- 分子量：200-500 Da
- LogP：-5至5
- H键捐助者：0-5
- H键受体：0-10

#### 药物规则()

```python
rule_of_drug(mol: Union[str, Chem.Mol]) -> bool
```

组合药物相似性标准。

**标准：**
- 通过五规则
- 通过 Veber 规则
- 无 PAINS 子结构

#### 黄金三角()

```python
golden_triangle(mol: Union[str, Chem.Mol]) -> bool
```

药物相似性平衡的金三角。

**标准：**
- 200≤MW≤50×LogP+400
- LogP：-2至5

#### pains_filter()

```python
pains_filter(mol: Union[str, Chem.Mol]) -> bool
```

泛测定干扰化合物 (PAINS) 过滤器。

**返回：** 如果分子不包含 PAINS 子结构，则返回 True

---

## 模块：medchem.structural

### 类：CommonAlertsFilters

筛选源自 ChEMBL 和文献的常见结构警报。

**构造函数：**
```python
CommonAlertsFilters()
```

**方法：**

```python
__call__(mols: List[Chem.Mol], n_jobs: int = 1, progress: bool = False) -> List[Dict]
```

将常见警报过滤器应用于分子列表。

**返回：** 带有键的字典列表：
- `has_alerts`：指示分子是否有警报的布尔值
- `alert_details`：匹配的警报模式列表
- `num_alerts`：找到的警报数量

```python
check_mol(mol: Chem.Mol) -> Tuple[bool, List[str]]
```

检查单个分子的结构警报。

**返回：** (has_alerts, list_of_alert_names) 元组

### 类：NIBR 过滤器

诺华 NIBR 药物化学过滤器。

**构造函数：**
```python
NIBRFilters()
```

**方法：**

```python
__call__(mols: List[Chem.Mol], n_jobs: int = 1, progress: bool = False) -> List[bool]
```

将 NIBR 过滤器应用于分子。

**返回：** 布尔值列表（如果分子通过则为 True）

### 类：LillyDemeritsFilters

礼来公司基于过失的结构性警报系统（275 条规则）。

**构造函数：**
```python
LillyDemeritsFilters()
```

**方法：**

```python
__call__(mols: List[Chem.Mol], n_jobs: int = 1, progress: bool = False) -> List[Dict]
```

计算分子的礼来缺点。

**返回：** 带有键的字典列表：
- `demerits`：总记过分
- `passes`：布尔值（如果缺点 ≤ 100，则为 True）
- `matched_patterns`：带分数的匹配模式列表

---

## 模块：medchem.function

用于常见操作的高级功能 API。

### nibr_filter()

```python
nibr_filter(mols: List[Chem.Mol], n_jobs: int = 1) -> List[bool]
```

使用功能 API 应用 NIBR 过滤器。

**参数：**
- `mols`：分子列表
- `n_jobs`：并行化级别

**返回：** 通过/失败布尔值列表

### common_alerts_filter()

```python
common_alerts_filter(mols: List[Chem.Mol], n_jobs: int = 1) -> List[Dict]
```

使用功能 API 应用常见警报过滤器。
**返回：**结果字典列表

### lilly_demerits_filter()

```python
lilly_demerits_filter(mols: List[Chem.Mol], n_jobs: int = 1) -> List[Dict]
```

使用功能 API 计算 Lilly 过失。

---

## 模块：medchem.groups

### 类别：化学组

检测分子中的特定化学基团。

**构造函数：**
```python
ChemicalGroup(groups: List[str], custom_smarts: Optional[Dict[str, str]] = None)
```

**参数：**
- `groups`：预定义组名称列表
- `custom_smarts`：将自定义组名称映射到 SMARTS 模式的字典

**预定义组：**
- `"hinge_binders"`：激酶铰链结合基序
- `"phosphate_binders"`：磷酸盐结合基团
- `"michael_acceptors"`：迈克尔受体亲电子试剂
- `"reactive_groups"`：一般反应功能

**方法：**

```python
has_match(mols: List[Chem.Mol]) -> List[bool]
```

检查分子是否包含任何指定的基团。

```python
get_matches(mol: Chem.Mol) -> Dict[str, List[Tuple]]
```

获取单个分子的详细匹配信息。

**返回：** 将组名称映射到原子索引列表的字典

```python
get_all_matches(mols: List[Chem.Mol]) -> List[Dict]
```

获取所有分子的匹配信息。

**示例：**
```python
group = mc.groups.ChemicalGroup(groups=["hinge_binders", "phosphate_binders"])
matches = group.get_all_matches(mol_list)
```

---

## 模块：medchem.catalogs

### 类：NamedCatalogs

访问精选的化学品目录。

**可用目录：**
- `"functional_groups"`：常见功能组
- `"protecting_groups"`：保护组结构
- `"reagents"`：常用试剂
- `"fragments"`：标准片段

**用途：**
```python
catalog = mc.catalogs.NamedCatalogs.get("functional_groups")
matches = catalog.get_matches(mol)
```

---

## 模块：medchem.complexity

计算分子复杂性指标。

### 计算复杂度()

```python
calculate_complexity(mol: Chem.Mol, method: str = "bertz") -> float
```

计算分子的复杂度分数。

**参数：**
- `mol`：RDKit 分子
- `method`：复杂度度量（“bertz”、“whitlock”、“barone”）

**返回：** 复杂度分数（越高=越复杂）

### 类：复杂性过滤器

按复杂度阈值过滤分子。

**构造函数：**
```python
ComplexityFilter(max_complexity: float, method: str = "bertz")
```

**方法：**

```python
__call__(mols: List[Chem.Mol], n_jobs: int = 1) -> List[bool]
```

过滤超过复杂性阈值的分子。

---

## 模块：medchem.constraints

### 类：约束

应用基于自定义属性的约束。

**构造函数：**
```python
Constraints(
    mw_range: Optional[Tuple[float, float]] = None,
    logp_range: Optional[Tuple[float, float]] = None,
    tpsa_max: Optional[float] = None,
    tpsa_range: Optional[Tuple[float, float]] = None,
    hbd_max: Optional[int] = None,
    hba_max: Optional[int] = None,
    rotatable_bonds_max: Optional[int] = None,
    rings_range: Optional[Tuple[int, int]] = None,
    aromatic_rings_max: Optional[int] = None,
)
```

**参数：** 所有参数都是可选的。仅指定所需的约束。

**方法：**

```python
__call__(mols: List[Chem.Mol], n_jobs: int = 1) -> List[Dict]
```

对分子施加约束。

**返回：** 带有键的字典列表：
- `passes`：指示是否所有约束都通过的布尔值
- `violations`：失败的约束名称列表

**示例：**
```python
constraints = mc.constraints.Constraints(
    mw_range=(200, 500),
    logp_range=(-2, 5),
    tpsa_max=140
)
results = constraints(mols=mol_list, n_jobs=-1)
```

---

## 模块：medchem.query

用于复杂过滤的查询语言。

### 解析()

```python
parse(query: str) -> Query
```

将 medchem 查询字符串解析为 Query 对象。

**查询语法：**
- 运算符：`AND`、`OR`、`NOT`
- 比较：`<`、`>`、`<=`、`>=`、`==`、`!=`
- 属性：`complexity`、`lilly_demerits`、`mw`、`logp`、`tpsa`
- 规则：`rule_of_five`、`rule_of_cns` 等。
- 过滤器：`common_alerts`、`nibr_filter`、`pains_filter`

**查询示例：**
```python
"rule_of_five AND NOT common_alerts"
"rule_of_cns AND complexity < 400"
"mw > 200 AND mw < 500 AND logp < 5"
"(rule_of_five OR rule_of_oprea) AND NOT pains_filter"
```

### 类：查询

**方法：**

```python
apply(mols: List[Chem.Mol], n_jobs: int = 1) -> List[bool]
```

将解析的查询应用于分子。

**示例：**
```python
query = mc.query.parse("rule_of_five AND NOT common_alerts")
results = query.apply(mols=mol_list, n_jobs=-1)
passing_mols = [mol for mol, passes in zip(mol_list, results) if passes]
```

---

## 模块：medchem.utils

用于处理分子的实用函数。

### 批处理()

```python
batch_process(
    mols: List[Chem.Mol],
    func: Callable,
    n_jobs: int = 1,
    progress: bool = False,
    batch_size: Optional[int] = None
) -> List
```

并行批次处理分子。

**参数：**
- `mols`：分子列表
- `func`：应用于每个分子的函数
- `n_jobs`：并行工作线程数
- `progress`：显示进度条
- `batch_size`：处理批次的大小

### standardize_mol()

```python
standardize_mol(mol: Chem.Mol) -> Chem.Mol
```

标准化分子表示（净化、中和电荷等）。

---

## 常见模式

### 模式：并行处理

所有过滤器都支持并行化：

```python
# Use all CPU cores
results = filter_object(mols=mol_list, n_jobs=-1, progress=True)

# Use specific number of cores
results = filter_object(mols=mol_list, n_jobs=4, progress=True)
```

### 模式：组合多个过滤器

```python
import medchem as mc

# Apply multiple filters
rule_filter = mc.rules.RuleFilters(rule_list=["rule_of_five"])
alert_filter = mc.structural.CommonAlertsFilters()
lilly_filter = mc.structural.LillyDemeritsFilters()

# Get results
rule_results = rule_filter(mols=mol_list, n_jobs=-1)
alert_results = alert_filter(mols=mol_list, n_jobs=-1)
lilly_results = lilly_filter(mols=mol_list, n_jobs=-1)

# Combine criteria
passing_mols = [
    mol for i, mol in enumerate(mol_list)
    if rule_results[i]["passes"]
    and not alert_results[i]["has_alerts"]
    and lilly_results[i]["passes"]
]
```

### 模式：使用 DataFrame

```python
import pandas as pd
import datamol as dm
import medchem as mc

# Load data
df = pd.read_csv("molecules.csv")
df["mol"] = df["smiles"].apply(dm.to_mol)

# Apply filters
rfilter = mc.rules.RuleFilters(rule_list=["rule_of_five", "rule_of_cns"])
results = rfilter(mols=df["mol"].tolist(), n_jobs=-1)

# Add results to dataframe
df["passes_ro5"] = [r["rule_of_five"] for r in results]
df["passes_cns"] = [r["rule_of_cns"] for r in results]

# Filter dataframe
filtered_df = df[df["passes_ro5"] & df["passes_cns"]]
```