<!-- 此文件由机器翻译自 smarts_patterns.md -->

# RDKit 的常见 SMARTS 模式

本文档提供了 RDKit 中用于子结构搜索的常用 SMARTS 模式的集合。

## 功能组

### 醇类

```python
# Primary alcohol
'[CH2][OH1]'

# Secondary alcohol
'[CH1]([OH1])[CH3,CH2]'

# Tertiary alcohol
'[C]([OH1])([C])([C])[C]'

# Any alcohol
'[OH1][C]'

# Phenol
'c[OH1]'
```

### 醛和酮

<<<代码块_1>>>

### 羧酸及其衍生物

<<<代码块_2>>>

### 胺类

<<<代码块_3>>>

### 醚类

<<<代码块_4>>>

### 卤化物

<<<代码块_5>>>

### 腈和硝基基团

<<<代码块_6>>>

### 硫醇和硫化物

```python
# Thiol
'[C][SH1]'

# Sulfide
'[C][S][C]'

# Disulfide
'[C][S][S][C]'

# Sulfoxide
'[C][S](=O)[C]'

# Sulfone
'[C][S](=O)(=O)[C]'
```

## 环形系统

### 简单戒指

```python
# Benzene ring
'c1ccccc1'
'[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'  # Explicit atoms

# Cyclohexane
'C1CCCCC1'

# Cyclopentane
'C1CCCC1'

# Any 3-membered ring
'[r3]'

# Any 4-membered ring
'[r4]'

# Any 5-membered ring
'[r5]'

# Any 6-membered ring
'[r6]'

# Any 7-membered ring
'[r7]'
```

### 芳香环

```python
# Aromatic carbon in ring
'[cR]'

# Aromatic nitrogen in ring (pyridine, etc.)
'[nR]'

# Aromatic oxygen in ring (furan, etc.)
'[oR]'

# Aromatic sulfur in ring (thiophene, etc.)
'[sR]'

# Any aromatic ring
'a1aaaaa1'
```

### 杂环

```python
# Pyridine
'n1ccccc1'

# Pyrrole
'n1cccc1'

# Furan
'o1cccc1'

# Thiophene
's1cccc1'

# Imidazole
'n1cncc1'

# Pyrimidine
'n1cnccc1'

# Thiazole
'n1ccsc1'

# Oxazole
'n1ccoc1'
```

### 熔接环

```python
# Naphthalene
'c1ccc2ccccc2c1'

# Indole
'c1ccc2[nH]ccc2c1'

# Quinoline
'n1cccc2ccccc12'

# Benzimidazole
'c1ccc2[nH]cnc2c1'

# Purine
'n1cnc2ncnc2c1'
```

### 大环化合物

```python
# Rings with 8 or more atoms
'[r{8-}]'

# Rings with 9-15 atoms
'[r{9-15}]'

# Rings with more than 12 atoms (macrocycles)
'[r{12-}]'
```

## 具体结构特点

### 脂肪族与芳香族

```python
# Aliphatic carbon
'[C]'

# Aromatic carbon
'[c]'

# Aliphatic carbon in ring
'[CR]'

# Aromatic carbon (alternative)
'[cR]'
```

### 立体化学

```python
# Tetrahedral center with clockwise chirality
'[C@]'

# Tetrahedral center with counterclockwise chirality
'[C@@]'

# Any chiral center
'[C@,C@@]'

# E double bond
'C/C=C/C'

# Z double bond
'C/C=C\\C'
```

### 杂交

```python
# SP hybridization (triple bond)
'[CX2]'

# SP2 hybridization (double bond or aromatic)
'[CX3]'

# SP3 hybridization (single bonds)
'[CX4]'
```

### 充电

```python
# Positive charge
'[+]'

# Negative charge
'[-]'

# Specific charge
'[+1]'
'[-1]'
'[+2]'

# Positively charged nitrogen
'[N+]'

# Negatively charged oxygen
'[O-]'

# Carboxylate anion
'C(=O)[O-]'

# Ammonium cation
'[N+]([C])([C])([C])[C]'
```

## 药效团特征

### 氢键供体

```python
# Hydroxyl
'[OH]'

# Amine
'[NH,NH2]'

# Amide NH
'[N][C](=O)'

# Any H-bond donor
'[OH,NH,NH2,NH3+]'
```

### 氢键受体

```python
# Carbonyl oxygen
'[O]=[C,S,P]'

# Ether oxygen
'[OX2]'

# Ester oxygen
'C(=O)[O]'

# Nitrogen acceptor
'[N;!H0]'

# Any H-bond acceptor
'[O,N]'
```

### 疏水基团

```python
# Alkyl chain (4+ carbons)
'CCCC'

# Branched alkyl
'C(C)(C)C'

# Aromatic rings (hydrophobic)
'c1ccccc1'
```

### 芳香相互作用

```python
# Benzene for pi-pi stacking
'c1ccccc1'

# Heterocycle for pi-pi
'[a]1[a][a][a][a][a]1'

# Any aromatic ring
'[aR]'
```

## 类似药物的片段

### 利宾斯基碎片

```python
# Aromatic ring with substituents
'c1cc(*)ccc1'

# Aliphatic chain
'CCCC'

# Ether linkage
'[C][O][C]'

# Amine (basic center)
'[N]([C])([C])'
```

### 常用脚手架

```python
# Benzamide
'c1ccccc1C(=O)N'

# Sulfonamide
'S(=O)(=O)N'

# Urea
'[N][C](=O)[N]'

# Guanidine
'[N]C(=[N])[N]'

# Phosphate
'P(=O)([O-])([O-])[O-]'
```

### 特权结构

```python
# Biphenyl
'c1ccccc1-c2ccccc2'

# Benzopyran
'c1ccc2OCCCc2c1'

# Piperazine
'N1CCNCC1'

# Piperidine
'N1CCCCC1'

# Morpholine
'N1CCOCC1'
```

## 反应基团

### 亲电子试剂

```python
# Acyl chloride
'C(=O)Cl'

# Alkyl halide
'[C][Cl,Br,I]'

# Epoxide
'C1OC1'

# Michael acceptor
'C=C[C](=O)'
```

### 亲核试剂

```python
# Primary amine
'[NH2][C]'

# Thiol
'[SH][C]'

# Alcohol
'[OH][C]'
```

## 毒性警报（疼痛）

```python
# Rhodanine
'S1C(=O)NC(=S)C1'

# Catechol
'c1ccc(O)c(O)c1'

# Quinone
'O=C1C=CC(=O)C=C1'

# Hydroquinone
'OC1=CC=C(O)C=C1'

# Alkyl halide (reactive)
'[C][I,Br]'

# Michael acceptor (reactive)
'C=CC(=O)[C,N]'
```

## 金属装订

```python
# Carboxylate (metal chelator)
'C(=O)[O-]'

# Hydroxamic acid
'C(=O)N[OH]'

# Catechol (iron chelator)
'c1c(O)c(O)ccc1'

# Thiol (metal binding)
'[SH]'

# Histidine-like (metal binding)
'c1ncnc1'
```

## 大小和复杂性过滤器

```python
# Long aliphatic chains (>6 carbons)
'CCCCCCC'

# Highly branched (quaternary carbon)
'C(C)(C)(C)C'

# Multiple rings
'[R]~[R]'  # Two rings connected

# Spiro center
'[C]12[C][C][C]1[C][C]2'
```

## 特殊图案

### 原子计数

```python
# Any atom
'[*]'

# Heavy atom (not H)
'[!H]'

# Carbon
'[C,c]'

# Heteroatom
'[!C;!H]'

# Halogen
'[F,Cl,Br,I]'
```

### 债券类型

```python
# Single bond
'C-C'

# Double bond
'C=C'

# Triple bond
'C#C'

# Aromatic bond
'c:c'

# Any bond
'C~C'
```

### 戒指会员资格

```python
# In any ring
'[R]'

# Not in ring
'[!R]'

# In exactly one ring
'[R1]'

# In exactly two rings
'[R2]'

# Ring bond
'[R]~[R]'
```

### 学位和连通性

```python
# Total degree 1 (terminal atom)
'[D1]'

# Total degree 2 (chain)
'[D2]'

# Total degree 3 (branch point)
'[D3]'

# Total degree 4 (highly branched)
'[D4]'

# Connected to exactly 2 carbons
'[C]([C])[C]'
```

## 用法示例

```python
from rdkit import Chem

# Create SMARTS query
pattern = Chem.MolFromSmarts('[CH2][OH1]')  # Primary alcohol

# Search molecule
mol = Chem.MolFromSmiles('CCO')
matches = mol.GetSubstructMatches(pattern)

# Multiple patterns
patterns = {
    'alcohol': '[OH1][C]',
    'amine': '[NH2,NH1][C]',
    'carboxylic_acid': 'C(=O)[OH1]'
}

# Check for functional groups
for name, smarts in patterns.items():
    query = Chem.MolFromSmarts(smarts)
    if mol.HasSubstructMatch(query):
        print(f"Found {name}")
```

## 编写 SMARTS 的技巧

1. **需要时具体化：** 使用原子属性 [CX3] 而不仅仅是 [C]
2. **为了清楚起见，使用括号：** [C] 与 C（芳香族）不同
3. **考虑芳香性：**小写字母（c、n、o）是芳香性的
4. **检查环成员资格：** [R] 表示环内，[!R] 表示非环内
5. **使用递归 SMARTS:** $(...) 用于复杂模式
6. **测试模式：** 始终在已知分子上验证 SMARTS
7. **从简单开始：** 逐步构建复杂的模式

## 常见 SMARTS 语法

- `[C]` - 脂肪族碳
- `[c]` - 芳香碳
- `[CX4]` - 具有 4 个连接的碳 (sp3)
- `[CX3]` - 具有 3 个连接的碳 (sp2)
- `[CX2]` - 具有 2 个连接的碳 (sp)
- `[CH3]` - 甲基
- `[R]` - 在环中
- `[r6]` - 在 6 元环中
- `[r{5-7}]` - 5、6 或 7 元环
- `[D2]` - 2 级（2 个邻居）
- `[+]` - 正电荷
- `[-]` - 负电荷
- `[!C]` - 非碳
- `[#6]` - 原子序数为 6 的元素（碳）
- `~` - 任何债券类型
- `-` - 单键
- `=` - 双键
- `#` - 三键
- `:` - 芳香键
- `@` - 顺时针手性
- `@@` - 逆时针手性