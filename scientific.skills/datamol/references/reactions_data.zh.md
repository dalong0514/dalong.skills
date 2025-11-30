<!-- 此文件由机器翻译自 reactions_data.md -->

# Datamol 反应和数据模块参考

## 反应模块 (`datamol.reactions`)

反应模块支持使用 SMARTS 反应模式以编程方式应用化学转化。

### 应用化学反应

#### `dm.reactions.apply_reaction(rxn, reactants, as_smiles=False, sanitize=True, single_product_group=True, rm_attach=True, product_index=0)`
对反应物分子进行化学反应。
- **参数**：
  - `rxn`：反应对象（来自 SMARTS 模式）
  - `reactants`：反应物分子元组
  - `as_smiles`：返回 SMILES 字符串（True）或分子对象（False）
  - `sanitize`：净化产品分子
  - `single_product_group`：返回单个产品 (True) 或所有产品组 (False)
  - `rm_attach`：删除附着点标记
  - `product_index`：从反应中返回哪个产品
- **返回**：产品分子或微笑
- **示例**：
  ```python
  from rdkit import Chem

  # Define reaction: alcohol + carboxylic acid → ester
  rxn = Chem.rdChemReactions.ReactionFromSmarts(
      '[C:1][OH:2].[C:3](=[O:4])[OH:5]>>[C:1][O:2][C:3](=[O:4])'
  )

  # Apply to reactants
  alcohol = dm.to_mol("CCO")
  acid = dm.to_mol("CC(=O)O")
  product = dm.reactions.apply_reaction(rxn, (alcohol, acid))
  ```

### 创建反应

反应通常是使用 RDKit 根据 SMARTS 模式创建的：
<<<代码块_1>>>

### 验证函数

该模块包含以下功能：
- **检查分子是否是反应物**：验证分子是否与反应物模式匹配
- **验证反应**：检查反应是否综合合理
- **处理反应文件**：从文件或数据库加载反应

### 常见反应模式

**酰胺形成**：
<<<代码块_2>>>

**铃木联轴器**：
<<<代码块_3>>>

**功能组转换**：
<<<代码块_4>>>

### 工作流程示例

<<<代码块_5>>>

### 关键概念

- **SMARTS**：SMiles ARbitrary 目标规范 - 反应模式语言
- **原子映射**：像 [C:1] 这样的数字通过反应保留原子身份
- **连接点**：[1*] 代表通用连接点
- **反应验证**：并非所有 SMARTS 反应在化学上都是合理的

---

## 数据模块 (`datamol.data`)

数据模块可以方便地访问精选的分子数据集以进行测试和学习。

### 可用数据集

#### `dm.data.cdk2(as_df=True, mol_column='mol')`
RDKit CDK2 数据集 - 激酶抑制剂数据。
- **参数**：
  - `as_df`：返回为 DataFrame (True) 或分子列表 (False)
  - `mol_column`：分子列的名称
- **返回**：具有分子结构和活性数据的数据集
- **用例**：用于算法测试的小数据集
- **示例**：
  <<<代码块_6>>>

#### `dm.data.freesolv()`
FreeSolv 数据集 - 实验和计算的水合自由能。
- **内容**：642 个分子，其中：
  - IUPAC 名称
  - 微笑字符串
  - 实验水合自由能值
  - 计算值
- **警告**：“仅用作教学和测试目的的玩具数据集”
- **不适合**：基准测试或生产模型培训
- **示例**：
  ```python
  freesolv_df = dm.data.freesolv()
  # Columns: iupac, smiles, expt (kcal/mol), calc (kcal/mol)
  ```

#### `dm.data.solubility(as_df=True, mol_column='mol')`
具有训练/测试分割的 RDKit 溶解度数据集。
- **内容**：具有预定义分割的水溶性数据
- **列**：包括带有“训练”或“测试”值的“拆分”列
- **用例**：通过正确的训练/测试分离来测试 ML 工作流程
- **示例**：
  ```python
  sol_df = dm.data.solubility(as_df=True)

  # Split into train/test
  train_df = sol_df[sol_df['split'] == 'train']
  test_df = sol_df[sol_df['split'] == 'test']

  # Use for model development
  X_train = dm.to_fp(train_df[mol_column])
  y_train = train_df['solubility']
  ```

### 使用指南

**用于测试和教程**：
```python
# Quick dataset for testing code
df = dm.data.cdk2()
mols = df['mol'].tolist()

# Test descriptor calculation
descriptors_df = dm.descriptors.batch_compute_many_descriptors(mols)

# Test clustering
clusters = dm.cluster_mols(mols, cutoff=0.3)
```

**对于学习工作流程**：
```python
# Complete ML pipeline example
sol_df = dm.data.solubility()

# Preprocessing
train = sol_df[sol_df['split'] == 'train']
test = sol_df[sol_df['split'] == 'test']

# Featurization
X_train = dm.to_fp(train['mol'])
X_test = dm.to_fp(test['mol'])

# Model training (example)
from sklearn.ensemble import RandomForestRegressor
model = RandomForestRegressor()
model.fit(X_train, train['solubility'])
predictions = model.predict(X_test)
```

### 重要提示

- **玩具数据集**：专为教学目的而设计，而非生产用途
- **小尺寸**：适合快速测试的化合物数量有限
- **预处理**：数据已清理和格式化
- **引用**：如果发布，请检查数据集文档以了解正确的归属

### 最佳实践

1. **仅用于开发**：不要从玩具数据集中得出科学结论
2. **验证真实数据**：始终在实际项目数据上测试生产代码
3. **正确归属**：如果在出版物中使用，请引用原始数据源
4. **了解限制**：了解每个数据集的范围和质量