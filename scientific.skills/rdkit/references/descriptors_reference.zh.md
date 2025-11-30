<!-- 此文件由机器翻译自 descriptors_reference.md -->

# RDKit 分子描述符参考

RDKit 的 `Descriptors` 模块中提供的分子描述符的完整参考。

## 用法

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

mol = Chem.MolFromSmiles('CCO')

# Calculate individual descriptor
mw = Descriptors.MolWt(mol)

# Calculate all descriptors at once
all_desc = Descriptors.CalcMolDescriptors(mol)
```

## 分子量和质量

### 摩尔重量
分子的平均分子量。
<<<代码块_1>>>

### 精确摩尔重量
使用同位素组成的精确分子量。
<<<代码块_2>>>

### 重原子摩尔重量
忽略氢的平均分子量。
<<<代码块_3>>>

## 亲脂性

### 分子日志
Wildman-Crippen LogP（辛醇-水分配系数）。
<<<代码块_4>>>

### 分子MR
Wildman-Crippen 摩尔折射率。
<<<代码块_5>>>

## 极地表面积

### TPSA
基于碎片贡献的拓扑极表面积 (TPSA)。
<<<代码块_6>>>

### 拉布特ASA
拉布特的近似表面积 (ASA)。
```python
Descriptors.LabuteASA(mol)
```

## 氢键

### NumHD 捐赠者
氢键供体的数量（N-H 和 O-H）。
```python
Descriptors.NumHDonors(mol)
```

### NumHAAcceptors
氢键受体（N 和 O）的数量。
```python
Descriptors.NumHAcceptors(mol)
```

### 没有计数
N 和 O 原子的数量。
```python
Descriptors.NOCount(mol)
```

### NHOH 计数
N-H 和 O-H 键的数量。
```python
Descriptors.NHOHCount(mol)
```

## 原子计数

### 重原子计数
重原子数（非氢）。
```python
Descriptors.HeavyAtomCount(mol)
```

### 杂原子数
杂原子数（非 C 和非 H）。
```python
Descriptors.NumHeteroatoms(mol)
```

### 价电子数
价电子总数。
```python
Descriptors.NumValenceElectrons(mol)
```

### NumRadicalElectrons
自由基电子数。
```python
Descriptors.NumRadicalElectrons(mol)
```

## 环描述符

### 环数
环数。
```python
Descriptors.RingCount(mol)
```

### NumAromaticRings
芳香环的数量。
```python
Descriptors.NumAromaticRings(mol)
```

### NumSaturatedRings
饱和环的数量。
```python
Descriptors.NumSaturatedRings(mol)
```

### NumAliphaticRings
脂肪族（非芳香族）环的数量。
```python
Descriptors.NumAliphaticRings(mol)
```

### NumAromaticCarbocycles
芳香碳环（仅含碳的环）的数量。
```python
Descriptors.NumAromaticCarbocycles(mol)
```

### NumAromaticHeterocycles
芳香杂环（含有杂原子的环）的数量。
```python
Descriptors.NumAromaticHeterocycles(mol)
```

### NumSaturatedCarbocycles
饱和碳环的数量。
```python
Descriptors.NumSaturatedCarbocycles(mol)
```

### NumSaturatedHeterocycles
饱和杂环的数量。
```python
Descriptors.NumSaturatedHeterocycles(mol)
```

### NumAliphaticCarbocycles
脂肪族碳环的数量。
```python
Descriptors.NumAliphaticCarbocycles(mol)
```

### NumAliphaticHeterocycles
脂肪族杂环的数量。
```python
Descriptors.NumAliphaticHeterocycles(mol)
```

## 可旋转债券

### NumRotableBonds
可旋转键的数量（灵活性）。
```python
Descriptors.NumRotatableBonds(mol)
```

## 芳香原子

### NumAromaticAtoms
芳香原子数。
```python
Descriptors.NumAromaticAtoms(mol)
```

## 分数描述符

### 分数Csp3
sp3杂化的碳分数。
```python
Descriptors.FractionCsp3(mol)
```

## 复杂度描述符

### 贝尔茨CT
Bertz 复杂度指数。
```python
Descriptors.BertzCT(mol)
```

### 工控机
信息内容（复杂性衡量）。
```python
Descriptors.Ipc(mol)
```

## Kappa 形状指数

基于图不变量的分子形状描述符。

###河童1
第一个 kappa 形状索引。
```python
Descriptors.Kappa1(mol)
```

###河童2
第二个 kappa 形状指数。
```python
Descriptors.Kappa2(mol)
```

###河童3
第三个 kappa 形状指数。
```python
Descriptors.Kappa3(mol)
```

## Chi 连接指数

分子连接指数。

### Chi0、Chi1、Chi2、Chi3、Chi4
简单的 chi 连接指数。
```python
Descriptors.Chi0(mol)
Descriptors.Chi1(mol)
Descriptors.Chi2(mol)
Descriptors.Chi3(mol)
Descriptors.Chi4(mol)
```

### Chi0n、Chi1n、Chi2n、Chi3n、Chi4n
价态修正的 chi 连接指数。
```python
Descriptors.Chi0n(mol)
Descriptors.Chi1n(mol)
Descriptors.Chi2n(mol)
Descriptors.Chi3n(mol)
Descriptors.Chi4n(mol)
```

### Chi0v、Chi1v、Chi2v、Chi3v、Chi4v
价chi连通性指数。
```python
Descriptors.Chi0v(mol)
Descriptors.Chi1v(mol)
Descriptors.Chi2v(mol)
Descriptors.Chi3v(mol)
Descriptors.Chi4v(mol)
```

## 霍尔-基尔阿尔法

### 霍尔基尔阿尔法
Hall-Kier α 值（分子灵活性）。
```python
Descriptors.HallKierAlpha(mol)
```

## 巴拉班 J 指数

### 巴拉班J
Balaban 的 J 索引（分支描述符）。
```python
Descriptors.BalabanJ(mol)
```

## EState 指数

电拓扑状态指数。

### 最大E状态索引
最大 E 状态值。
```python
Descriptors.MaxEStateIndex(mol)
```

### MinEStateIndex
最小 E 状态值。
```python
Descriptors.MinEStateIndex(mol)
```

### 最大AbsEStateIndex
最大绝对 E 状态值。
```python
Descriptors.MaxAbsEStateIndex(mol)
```

### 最小AbsEStateIndex
最小绝对 E 状态值。
```python
Descriptors.MinAbsEStateIndex(mol)
```

## 部分收费

### 最大部分充电
最大部分电荷。
```python
Descriptors.MaxPartialCharge(mol)
```

### 最小部分充电
最低部分费用。
```python
Descriptors.MinPartialCharge(mol)
```

### 最大绝对部分电荷
最大绝对部分电荷。
```python
Descriptors.MaxAbsPartialCharge(mol)
```

### MinAbsPartialCharge
最小绝对部分电荷。
```python
Descriptors.MinAbsPartialCharge(mol)
```

## 指纹密度

测量分子指纹的密度。

### FpDensityMorgan1
半径 1 处的摩根指纹密度。
```python
Descriptors.FpDensityMorgan1(mol)
```

### FpDensityMorgan2
半径 2 处的摩根指纹密度。
```python
Descriptors.FpDensityMorgan2(mol)
```

### FpDensityMorgan3
半径 3 处的摩根指纹密度。
```python
Descriptors.FpDensityMorgan3(mol)
```
## PEOE VSA 描述符

轨道电负性部分均衡 (PEOE) VSA 描述符。

### PEOE_VSA1 到 PEOE_VSA14
使用部分电荷和表面积贡献的 MOE 类型描述符。
```python
Descriptors.PEOE_VSA1(mol)
# ... through PEOE_VSA14
```

## SMR VSA 描述符

分子折射率 VSA 描述符。

### SMR_VSA1 到 SMR_VSA10
使用 MR 贡献和表面积的 MOE 类型描述符。
```python
Descriptors.SMR_VSA1(mol)
# ... through SMR_VSA10
```

## SLogP VSA 描述符

LogP VSA 描述符。

### SLogP_VSA1 到 SLogP_VSA12
使用 LogP 贡献和表面积的 MOE 类型描述符。
```python
Descriptors.SLogP_VSA1(mol)
# ... through SLogP_VSA12
```

## EState VSA 描述符

### EState_VSA1 到 EState_VSA11
使用 E 状态指数和表面积的 MOE 类型描述符。
```python
Descriptors.EState_VSA1(mol)
# ... through EState_VSA11
```

## VSA 描述符

范德华表面积描述符。

### VSA_EState1 到 VSA_EState10
EState VSA 描述符。
```python
Descriptors.VSA_EState1(mol)
# ... through VSA_EState10
```

## BCUT 描述符

Burden-CAS-德克萨斯大学特征值描述符。

### BCUT2D_MWHI
按分子量加权的负荷矩阵的最高特征值。
```python
Descriptors.BCUT2D_MWHI(mol)
```

### BCUT2D_MWLOW
按分子量加权的负荷矩阵的最低特征值。
```python
Descriptors.BCUT2D_MWLOW(mol)
```

### BCUT2D_CHGHI
按部分电荷加权的最高特征值。
```python
Descriptors.BCUT2D_CHGHI(mol)
```

### BCUT2D_CHGLO
按部分电荷加权的最低特征值。
```python
Descriptors.BCUT2D_CHGLO(mol)
```

### BCUT2D_LOGPHI
由 LogP 加权的最高特征值。
```python
Descriptors.BCUT2D_LOGPHI(mol)
```

### BCUT2D_LOGPLOW
由 LogP 加权的最低特征值。
```python
Descriptors.BCUT2D_LOGPLOW(mol)
```

### BCUT2D_MRHI
按摩尔折射率加权的最高特征值。
```python
Descriptors.BCUT2D_MRHI(mol)
```

### BCUT2D_MRLOW
按摩尔折射率加权的最低特征值。
```python
Descriptors.BCUT2D_MRLOW(mol)
```

## 自相关描述符

### AUTOCORR2D
2D 自相关描述符（如果启用）。
测量属性空间分布的各种自相关指数。

## MQN 描述符

分子量子数 - 42 个简单描述符。

### mqn1 到 mqn42
计算各种分子特征的整数描述符。
```python
# Access via CalcMolDescriptors
desc = Descriptors.CalcMolDescriptors(mol)
mqns = {k: v for k, v in desc.items() if k.startswith('mqn')}
```

## 量子电动力学

### qed
药物相似性的定量估计。
```python
Descriptors.qed(mol)
```

## 利宾斯基五法则

使用 Lipinski 的标准检查药物相似性：

```python
def lipinski_rule_of_five(mol):
    mw = Descriptors.MolWt(mol) <= 500
    logp = Descriptors.MolLogP(mol) <= 5
    hbd = Descriptors.NumHDonors(mol) <= 5
    hba = Descriptors.NumHAcceptors(mol) <= 10
    return mw and logp and hbd and hba
```

## 批量描述符计算

一次性计算所有描述符：

```python
from rdkit import Chem
from rdkit.Chem import Descriptors

mol = Chem.MolFromSmiles('CCO')

# Get all descriptors as dictionary
all_descriptors = Descriptors.CalcMolDescriptors(mol)

# Access specific descriptor
mw = all_descriptors['MolWt']
logp = all_descriptors['MolLogP']

# Get list of available descriptor names
from rdkit.Chem import Descriptors
descriptor_names = [desc[0] for desc in Descriptors._descList]
```

## 描述符类别摘要

1. **物理化学**：MolWt、MolLogP、MolMR、TPSA
2. **拓扑**：BertzCT、BalabanJ、Kappa 指数
3. **电子**：部分收费、电子状态指数
4. **形状**：Kappa 指数、BCUT 描述符
5. **连通性**：Chi 指数
6. **2D 指纹**：FpDensity 描述符
7. **原子计数**：重原子、杂原子、环
8. **药物相似性**：QED、Lipinski 参数
9. **灵活性**：NumRotatableBonds、HallKierAlpha
10. **表面积**：基于 VSA 的描述符

## 常见用例

### 药物相似性筛选

```python
def screen_druglikeness(mol):
    return {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'HBD': Descriptors.NumHDonors(mol),
        'HBA': Descriptors.NumHAcceptors(mol),
        'TPSA': Descriptors.TPSA(mol),
        'RotBonds': Descriptors.NumRotatableBonds(mol),
        'AromaticRings': Descriptors.NumAromaticRings(mol),
        'QED': Descriptors.qed(mol)
    }
```

### 类铅过滤

```python
def is_leadlike(mol):
    mw = 250 <= Descriptors.MolWt(mol) <= 350
    logp = Descriptors.MolLogP(mol) <= 3.5
    rot_bonds = Descriptors.NumRotatableBonds(mol) <= 7
    return mw and logp and rot_bonds
```

### 多样性分析

```python
def molecular_complexity(mol):
    return {
        'BertzCT': Descriptors.BertzCT(mol),
        'NumRings': Descriptors.RingCount(mol),
        'NumRotBonds': Descriptors.NumRotatableBonds(mol),
        'FractionCsp3': Descriptors.FractionCsp3(mol),
        'NumAromaticRings': Descriptors.NumAromaticRings(mol)
    }
```

## 提示

1. **对多个描述符使用批量计算**，避免冗余计算
2. **检查 None** - 某些描述符可能会为无效分子返回 None
3. **机器学习应用的描述符标准化**
4. **选择相关描述符** - 并非所有 200 多个描述符都对每个任务都有用
5. **单独考虑3D描述符**（需要3D坐标）
6. **验证范围** - 检查描述符值是否在预期范围内