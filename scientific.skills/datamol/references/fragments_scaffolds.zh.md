<!-- 此文件由机器翻译自 fragments_scaffolds.md -->

# Datamol 片段和支架参考

## 脚手架模块 (`datamol.scaffold`)

支架代表分子的核心结构，可用于识别结构族和分析构效关系 (SAR)。

### Murcko 脚手架

#### `dm.to_scaffold_murcko(mol)`
提取 Bemis-Murcko 支架（分子框架）。
- **方法**：去除侧链、固定环系统和连接体
- **返回**：代表支架的分子对象
- **用例**：识别化合物系列的核心结构
- **示例**：
  ```python
  mol = dm.to_mol("c1ccc(cc1)CCN")  # Phenethylamine
  scaffold = dm.to_scaffold_murcko(mol)
  scaffold_smiles = dm.to_smiles(scaffold)
  # Returns: 'c1ccccc1CC' (benzene ring + ethyl linker)
  ```

**支架分析工作流程**：
<<<代码块_1>>>

### 模糊支架

#### `dm.scaffold.fuzzy_scaffolding(mol, ...)`
生成具有必须出现在核心中的可执行组的模糊支架。
- **目的**：更灵活的支架定义，允许指定的功能组
- **用例**：Murcko 规则之外的自定义脚手架定义

### 应用程序

**基于支架的分割**（用于 ML 模型验证）：
<<<代码块_2>>>

**SAR分析**：
<<<代码块_3>>>

---

## 片段模块 (`datamol.fragment`)

分子破碎根据化学规则将分子分解成更小的碎片，这对于基于片段的药物设计和子结构分析非常有用。

### 金砖国家分裂

#### `dm.fragment.brics(mol, ...)`
使用 BRICS（打破逆合成有趣的化学子结构）片段分子。
- **方法**：根据 16 种具有化学意义的键类型进行剖析
- **考虑**：考虑化学环境和周围的子结构
- **返回**：片段 SMILES 字符串集
- **用例**：逆合成分析，基于片段的设计
- **示例**：
  <<<代码块_4>>>

### 回顾碎片

#### `dm.fragment.recap(mol, ...)`
使用 RECAP（逆合成组合分析程序）进行分子片段化。
- **方法**：根据 11 种预定义的键类型进行剖析
- **规则**：
  - 保留小于 5 个碳的烷基完整
  - 保留循环债券
- **返回**：片段 SMILES 字符串集
- **用例**：组合库设计
- **示例**：
  <<<代码块_5>>>

### MMPA 碎片

#### `dm.fragment.mmpa_frag(mol, ...)`
用于匹配分子对分析的片段。
- **目的**：生成适合识别分子对的片段
- **用例**：分析微小的结构变化如何影响属性
- **示例**：
  <<<代码块_6>>>

### 方法比较

|方法|债券类型 |保存周期|最适合 |
|--------|---------|------------------|----------|
|金砖国家 | 16 | 16是的 |逆合成分析、片段重组 |
|回顾 | 11 | 11是的 |组合文库设计 |
| MMPA |变量|取决于 |构效关系分析 |

### 碎片工作流程

```python
import datamol as dm

# 1. Fragment a molecule
mol = dm.to_mol("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
brics_frags = dm.fragment.brics(mol)
recap_frags = dm.fragment.recap(mol)

# 2. Analyze fragment frequency across library
all_fragments = []
for mol in molecule_library:
    frags = dm.fragment.brics(mol)
    all_fragments.extend(frags)

# 3. Identify common fragments
from collections import Counter
fragment_counts = Counter(all_fragments)
common_fragments = fragment_counts.most_common(20)

# 4. Convert fragments back to molecules (remove attachment points)
def clean_fragment(frag_smiles):
    # Remove [1*], [2*], etc. attachment point markers
    clean = frag_smiles.replace('[1*]', '[H]')
    return dm.to_mol(clean)
```

### 高级：基于片段的虚拟筛选

```python
# Build fragment library from known actives
active_fragments = set()
for active_mol in active_compounds:
    frags = dm.fragment.brics(active_mol)
    active_fragments.update(frags)

# Screen compounds for presence of active fragments
def score_by_fragments(mol, fragment_set):
    mol_frags = dm.fragment.brics(mol)
    overlap = mol_frags.intersection(fragment_set)
    return len(overlap) / len(mol_frags)

# Score screening library
scores = [score_by_fragments(mol, active_fragments) for mol in screening_lib]
```

### 关键概念

- **附着点**：在 SMILES 片段中标有 [1*]、[2*] 等
- **逆合成**：碎片模拟合成断开
- **具有化学意义**：典型合成键发生断裂
- **重组**：理论上片段可以重组成有效的分子