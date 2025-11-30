<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 医学化学
描述：“药物化学过滤器。应用药物相似性规则（Lipinski、Veber）、PAINS 过滤器、结构警报、复杂性指标，以进行化合物优先级排序和库过滤。”
---

# 医学化学

## 概述

Medchem 是一个用于药物发现工作流程中的分子过滤和优先级排序的 Python 库。应用数百个完善的新型分子过滤器、结构警报和药物化学规则，大规模有效地对化合物库进行分类和优先排序。规则和过滤器是特定于上下文的 - 与领域专业知识相结合作为指南。

## 何时使用此技能

该技能应该在以下情况下使用：
- 将药物相似规则（Lipinski、Veber 等）应用于化合物库
- 通过结构警报或 PAINS 模式过滤分子
- 确定先导化合物优化的优先顺序
- 评估化合物质量和药物化学特性
- 检测反应性或有问题的官能团
- 计算分子复杂性指标

## 安装

```bash
uv pip install medchem
```

## 核心能力

### 1. 药物化学规则

使用 `medchem.rules` 模块将已建立的药物相似性规则应用于分子。

**可用规则：**
- 五法则（利平斯基）
- 歌剧规则
- CNS规则
- 铅状规则（软和严格）
- 三法则
- Reos规则
- 药物规则
- 韦伯规则
- 金三角
- 痛苦过滤器

**单一规则应用：**

<<<代码块_1>>>

**带有规则过滤器的多个规则：**

<<<代码块_2>>>

**结果格式：**
结果以字典形式返回，其中包含通过/失败状态以及每条规则的详细信息。

### 2. 结构警报过滤器

使用 `medchem.structural` 模块检测潜在有问题的结构模式。

**可用过滤器：**

1. **常见警报** - 来自 ChEMBL 管理和文献的一般结构警报
2. **NIBR 过滤器** - 诺华生物医学研究所过滤器套件
3. **Lilly Demerits** - Eli Lilly 的基于过失的系统（275 条规则，超过 100 个过失的分子被拒绝）

**常见警报：**

<<<代码块_3>>>

**NIBR 过滤器：**

<<<代码块_4>>>

**莉莉的缺点：**

<<<代码块_5>>>

### 3. 用于高级操作的函数式 API

`medchem.functional` 模块为常见工作流程提供了便捷的功能。

**快速过滤：**

<<<代码块_6>>>

### 4. 化学基团检测

使用`medchem.groups`识别特定的化学基团和官能团。

**可用组：**
- 铰链活页夹
- 磷酸盐结合剂
- 迈克尔接受者
- 反应基团
- 自定义 SMARTS 模式

**用途：**

```python
import medchem as mc

# Create group detector
group = mc.groups.ChemicalGroup(groups=["hinge_binders"])

# Check for matches
has_matches = group.has_match(mol_list)

# Get detailed match information
matches = group.get_matches(mol)
```

### 5. 命名目录

通过 `medchem.catalogs` 访问精选的化学结构集合。

**可用目录：**
- 功能组
- 保护群体
- 常用试剂
- 标准片段

**用途：**

```python
import medchem as mc

# Access named catalogs
catalogs = mc.catalogs.NamedCatalogs

# Use catalog for matching
catalog = catalogs.get("functional_groups")
matches = catalog.get_matches(mol)
```

### 6.分子复杂性

使用 `medchem.complexity` 计算近似综合可访问性的复杂性指标。

**通用指标：**
- 贝尔茨复杂度
- 惠特洛克复杂性
- 男爵复杂性

**用途：**

```python
import medchem as mc

# Calculate complexity
complexity_score = mc.complexity.calculate_complexity(mol)

# Filter by complexity threshold
complex_filter = mc.complexity.ComplexityFilter(max_complexity=500)
results = complex_filter(mols=mol_list)
```

### 7. 约束过滤

使用 `medchem.constraints` 应用基于自定义属性的约束。

**约束示例：**
- 分子量范围
- LogP 界限
- TPSA 限制
- 可旋转债券数量

**用途：**

```python
import medchem as mc

# Define constraints
constraints = mc.constraints.Constraints(
    mw_range=(200, 500),
    logp_range=(-2, 5),
    tpsa_max=140,
    rotatable_bonds_max=10
)

# Apply constraints
results = constraints(mols=mol_list, n_jobs=-1)
```

### 8. Medchem 查询语言

使用专门的查询语言来实现复杂的过滤条件。

**查询示例：**
```
# Molecules passing Ro5 AND not having common alerts
"rule_of_five AND NOT common_alerts"

# CNS-like molecules with low complexity
"rule_of_cns AND complexity < 400"

# Leadlike molecules without Lilly demerits
"rule_of_leadlike AND lilly_demerits == 0"
```

**用途：**

```python
import medchem as mc

# Parse and apply query
query = mc.query.parse("rule_of_five AND NOT common_alerts")
results = query.apply(mols=mol_list, n_jobs=-1)
```

## 工作流程模式

### 模式 1：化合物库的初始分类

过滤大量化合物集合以识别类药物候选物。

```python
import datamol as dm
import medchem as mc
import pandas as pd

# Load compound library
df = pd.read_csv("compounds.csv")
mols = [dm.to_mol(smi) for smi in df["smiles"]]

# Apply primary filters
rule_filter = mc.rules.RuleFilters(rule_list=["rule_of_five", "rule_of_veber"])
rule_results = rule_filter(mols=mols, n_jobs=-1, progress=True)

# Apply structural alerts
alert_filter = mc.structural.CommonAlertsFilters()
alert_results = alert_filter(mols=mols, n_jobs=-1, progress=True)

# Combine results
df["passes_rules"] = rule_results["pass"]
df["has_alerts"] = alert_results["has_alerts"]
df["drug_like"] = df["passes_rules"] & ~df["has_alerts"]

# Save filtered compounds
filtered_df = df[df["drug_like"]]
filtered_df.to_csv("filtered_compounds.csv", index=False)
```

### 模式 2：先导优化过滤

在先导化合物优化过程中应用更严格的标准。

```python
import medchem as mc

# Create comprehensive filter
filters = {
    "rules": mc.rules.RuleFilters(rule_list=["rule_of_leadlike_strict"]),
    "alerts": mc.structural.NIBRFilters(),
    "lilly": mc.structural.LillyDemeritsFilters(),
    "complexity": mc.complexity.ComplexityFilter(max_complexity=400)
}

# Apply all filters
results = {}
for name, filt in filters.items():
    results[name] = filt(mols=candidate_mols, n_jobs=-1)

# Identify compounds passing all filters
passes_all = all(r["pass"] for r in results.values())
```

### 模式 3：识别特定化学基团

寻找含有特定官能团或支架的分子。

```python
import medchem as mc

# Create group detector for multiple groups
group_detector = mc.groups.ChemicalGroup(
    groups=["hinge_binders", "phosphate_binders"]
)

# Screen library
matches = group_detector.get_all_matches(mol_list)

# Filter molecules with desired groups
mol_with_groups = [mol for mol, match in zip(mol_list, matches) if match]
```

## 最佳实践

1. **上下文很重要**：不要盲目应用过滤器。了解生物靶标和化学空间。

2. **组合多个过滤器**：将规则、结构警报和领域知识结合使用以做出更好的决策。

3. **使用并行化**：对于大型数据集（>1000 个分子），始终使用 `n_jobs=-1` 进行并行处理。

4. **迭代细化**：从广泛的过滤器 (Ro5) 开始，然后根据需要应用更具体的标准（CNS、铅状）。
5. **文档过滤决策**：跟踪哪些分子被过滤掉以及过滤的原因以实现可重复性。

6. **验证结果**：请记住，上市药物通常无法通过标准过滤器 - 将这些作为指导方针，而不是绝对规则。

7. **考虑前药**：设计为前药的分子可能故意违反标准药物化学规则。

## 资源

### 参考资料/api_guide.md
全面的 API 参考，涵盖所有 medchem 模块，包含详细的函数签名、参数和返回类型。

### 参考文献/rules_catalog.md
可用规则、过滤器和警报的完整目录，以及描述、阈值和文献参考。

### 脚本/filter_molecules.py
用于批量过滤工作流程的生产就绪脚本。支持多种输入格式（CSV、SDF、SMILES）、可配置的过滤器组合和详细的报告。

**用途：**
```bash
python scripts/filter_molecules.py input.csv --rules rule_of_five,rule_of_cns --alerts nibr --output filtered.csv
```

## 文档

官方文档：https://medchem-docs.datamol.io/
GitHub 存储库：https://github.com/datamol-io/medchem