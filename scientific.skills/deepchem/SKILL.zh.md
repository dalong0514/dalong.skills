<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：深化学
描述：“分子机器学习工具包。属性预测（ADMET、毒性）、GNN（GCN、MPNN）、MoleculeNet 基准、预训练模型、特征化、用于药物发现 ML。”
---

# 深化学

## 概述

DeepChem 是一个综合性 Python 库，用于将机器学习应用于化学、材料科学和生物学。通过专门的神经网络、分子特征化方法和预训练模型，实现分子特性预测、药物发现、材料设计和生物分子分析。

## 何时使用此技能

该技能应该在以下情况下使用：
- 加载和处理分子数据（SMILES 字符串、SDF 文件、蛋白质序列）
- 预测分子特性（溶解度、毒性、结合亲和力、ADMET 特性）
- 化学/生物数据集的训练模型
- 使用 MoleculeNet 基准数据集（Tox21、BBBP、Delaney 等）
- 将分子转换为 ML 就绪特征（指纹、图形表示、描述符）
- 实现分子图神经网络（GCN、GAT、MPNN、AttentiveFP）
- 通过预训练模型应用迁移学习（ChemBERTa、GROVER、MolFormer）
- 预测晶体/材料特性（带隙、形成能）
- 分析蛋白质或DNA序列

## 核心能力

### 1. 分子数据加载和处理

DeepChem 为各种化学数据格式提供专门的加载器：

```python
import deepchem as dc

# Load CSV with SMILES
featurizer = dc.feat.CircularFingerprint(radius=2, size=2048)
loader = dc.data.CSVLoader(
    tasks=['solubility', 'toxicity'],
    feature_field='smiles',
    featurizer=featurizer
)
dataset = loader.create_dataset('molecules.csv')

# Load SDF files
loader = dc.data.SDFLoader(tasks=['activity'], featurizer=featurizer)
dataset = loader.create_dataset('compounds.sdf')

# Load protein sequences
loader = dc.data.FASTALoader()
dataset = loader.create_dataset('proteins.fasta')
```

**密钥加载器**：
- `CSVLoader`：带有分子标识符的表格数据
- `SDFLoader`：分子结构文件
- `FASTALoader`：蛋白质/DNA 序列
- `ImageLoader`：分子图像
- `JsonLoader`：JSON 格式的数据集

### 2. 分子特征化

将分子转换为 ML 模型的数字表示。

#### 特征器选择的决策树

<<<代码块_1>>>

#### 特征化示例

<<<代码块_2>>>

**选型指南**：
- **小数据集 (<1K)**：CircularFingerprint 或 RDKitDescriptors
- **中型数据集 (1K-100K)**：CircularFingerprint 或图形特征器
- **大型数据集 (>100K)**：图特征器（MolGraphConvFeaturizer、DMPNNFeaturizer）
- **迁移学习**：预训练模型特征化器 (GroverFeaturizer)

请参阅 `references/api_reference.md` 以获取完整的特征化器文档。

### 3.数据分割

**关键**：对于药物发现任务，请使用 `ScaffoldSplitter` 来防止训练和测试集中出现的相似分子结构的数据泄漏。

<<<代码块_3>>>

**可用的分离器**：
- `ScaffoldSplitter`：通过分子支架分割（防止泄漏）
- `ButinaSplitter`：基于聚类的分子分裂
- `MaxMinSplitter`：最大化集合之间的多样性
- `RandomSplitter`：随机分割
- `RandomStratifiedSplitter`：保留类分布

### 4.模型选择和训练

#### 快速选型指南

|数据集大小 |任务|推荐型号 |特征化器 |
|------------|---------|--------------------|------------|
| < 1K 样本 |任何 | SklearnModel（随机森林）|圆形指纹|
| 1K-100K |分类/回归| GBDT模型或多任务回归器|圆形指纹|
| > 10 万 |分子性质| GCN模型、AttentiveFP模型、DMPNN模型 | MolGraphConvFeaturizer | MolGraphConvFeaturizer |
|任意（较小的优先）|迁移学习 | ChemBERTa、GROVER、MolFormer |特定型号|
|晶体结构|材料特性| CGCNN模型、MEGNet模型 |基于结构|
|蛋白质序列|蛋白质特性 |普特伯特 |基于序列|

#### 示例：传统机器学习
<<<代码块_4>>>

#### 示例：深度学习
<<<代码块_5>>>

#### 示例：图神经网络
<<<代码块_6>>>

### 5. MoleculeNet 基准测试

快速访问 30 多个精选的基准数据集以及标准化的训练/有效/测试拆分：

```python
# Load benchmark dataset
tasks, datasets, transformers = dc.molnet.load_tox21(
    featurizer='GraphConv',  # or 'ECFP', 'Weave', 'Raw'
    splitter='scaffold',     # or 'random', 'stratified'
    reload=False
)
train, valid, test = datasets

# Train and evaluate
model = dc.models.GCNModel(n_tasks=len(tasks), mode='classification')
model.fit(train, nb_epoch=50)

metric = dc.metrics.Metric(dc.metrics.roc_auc_score)
test_score = model.evaluate(test, [metric])
```

**通用数据集**：
- **分类**：`load_tox21()`、`load_bbbp()`、`load_hiv()`、`load_clintox()`
- **回归**：`load_delaney()`、`load_freesolv()`、`load_lipo()`
- **量子属性**：`load_qm7()`、`load_qm8()`、`load_qm9()`
- **材料**：`load_perovskite()`、`load_bandgap()`、`load_mp_formation_energy()`

请参阅 `references/api_reference.md` 以获取完整的数据集列表。

### 6. 迁移学习

利用预训练模型来提高性能，尤其是在小型数据集上：

```python
# ChemBERTa (BERT pretrained on 77M molecules)
model = dc.models.HuggingFaceModel(
    model='seyonec/ChemBERTa-zinc-base-v1',
    task='classification',
    n_tasks=1,
    learning_rate=2e-5  # Lower LR for fine-tuning
)
model.fit(train, nb_epoch=10)

# GROVER (graph transformer pretrained on 10M molecules)
model = dc.models.GroverModel(
    task='regression',
    n_tasks=1
)
model.fit(train, nb_epoch=20)
```

**何时使用迁移学习**：
- 小数据集（< 1000 个样本）
- 新型分子支架
- 计算资源有限
- 需要快速原型制作

使用 `scripts/transfer_learning.py` 脚本来引导迁移学习工作流程。

### 7. 模型评估

```python
# Define metrics
classification_metrics = [
    dc.metrics.Metric(dc.metrics.roc_auc_score, name='ROC-AUC'),
    dc.metrics.Metric(dc.metrics.accuracy_score, name='Accuracy'),
    dc.metrics.Metric(dc.metrics.f1_score, name='F1')
]

regression_metrics = [
    dc.metrics.Metric(dc.metrics.r2_score, name='R²'),
    dc.metrics.Metric(dc.metrics.mean_absolute_error, name='MAE'),
    dc.metrics.Metric(dc.metrics.root_mean_squared_error, name='RMSE')
]

# Evaluate
train_scores = model.evaluate(train, classification_metrics)
test_scores = model.evaluate(test, classification_metrics)
```

### 8. 做出预测

```python
# Predict on test set
predictions = model.predict(test)

# Predict on new molecules
new_smiles = ['CCO', 'c1ccccc1', 'CC(C)O']
new_features = featurizer.featurize(new_smiles)
new_dataset = dc.data.NumpyDataset(X=new_features)

# Apply same transformations as training
for transformer in transformers:
    new_dataset = transformer.transform(new_dataset)

predictions = model.predict(new_dataset)
```

## 典型工作流程

### 工作流程 A：快速基准评估

为了在标准基准上评估模型：

```python
import deepchem as dc

# 1. Load benchmark
tasks, datasets, _ = dc.molnet.load_bbbp(
    featurizer='GraphConv',
    splitter='scaffold'
)
train, valid, test = datasets

# 2. Train model
model = dc.models.GCNModel(n_tasks=len(tasks), mode='classification')
model.fit(train, nb_epoch=50)

# 3. Evaluate
metric = dc.metrics.Metric(dc.metrics.roc_auc_score)
test_score = model.evaluate(test, [metric])
print(f"Test ROC-AUC: {test_score}")
```

### 工作流程 B：自定义数据预测

对于自定义分子数据集的训练：

```python
import deepchem as dc

# 1. Load and featurize data
featurizer = dc.feat.CircularFingerprint(radius=2, size=2048)
loader = dc.data.CSVLoader(
    tasks=['activity'],
    feature_field='smiles',
    featurizer=featurizer
)
dataset = loader.create_dataset('my_molecules.csv')

# 2. Split data (use ScaffoldSplitter for molecules!)
splitter = dc.splits.ScaffoldSplitter()
train, valid, test = splitter.train_valid_test_split(dataset)

# 3. Normalize (optional but recommended)
transformers = [dc.trans.NormalizationTransformer(
    transform_y=True, dataset=train
)]
for transformer in transformers:
    train = transformer.transform(train)
    valid = transformer.transform(valid)
    test = transformer.transform(test)

# 4. Train model
model = dc.models.MultitaskRegressor(
    n_tasks=1,
    n_features=2048,
    layer_sizes=[1000, 500],
    dropouts=0.25
)
model.fit(train, nb_epoch=50)

# 5. Evaluate
metric = dc.metrics.Metric(dc.metrics.r2_score)
test_score = model.evaluate(test, [metric])
```

### 工作流程 C：小数据集上的迁移学习

为了利用预训练模型：

```python
import deepchem as dc

# 1. Load data (pretrained models often need raw SMILES)
loader = dc.data.CSVLoader(
    tasks=['activity'],
    feature_field='smiles',
    featurizer=dc.feat.DummyFeaturizer()  # Model handles featurization
)
dataset = loader.create_dataset('small_dataset.csv')

# 2. Split data
splitter = dc.splits.ScaffoldSplitter()
train, test = splitter.train_test_split(dataset)

# 3. Load pretrained model
model = dc.models.HuggingFaceModel(
    model='seyonec/ChemBERTa-zinc-base-v1',
    task='classification',
    n_tasks=1,
    learning_rate=2e-5
)

# 4. Fine-tune
model.fit(train, nb_epoch=10)

# 5. Evaluate
predictions = model.predict(test)
```

请参阅 `references/workflows.md` 了解 8 个详细的工作流程示例，涵盖分子生成、材料科学、蛋白质分析等。

## 示例脚本

此技能包括 `scripts/` 目录中的三个生产就绪脚本：

### 1.`predict_solubility.py`
训练和评估溶解度预测模型。适用于 Delaney 基准或自定义 CSV 数据。

```bash
# Use Delaney benchmark
python scripts/predict_solubility.py

# Use custom data
python scripts/predict_solubility.py \
    --data my_data.csv \
    --smiles-col smiles \
    --target-col solubility \
    --predict "CCO" "c1ccccc1"
```

### 2.`graph_neural_network.py`
在分子数据上训练各种图神经网络架构。

```bash
# Train GCN on Tox21
python scripts/graph_neural_network.py --model gcn --dataset tox21

# Train AttentiveFP on custom data
python scripts/graph_neural_network.py \
    --model attentivefp \
    --data molecules.csv \
    --task-type regression \
    --targets activity \
    --epochs 100
```

### 3.`transfer_learning.py`
针对分子特性预测任务微调预训练模型（ChemBERTa、GROVER）。

```bash
# Fine-tune ChemBERTa on BBBP
python scripts/transfer_learning.py --model chemberta --dataset bbbp

# Fine-tune GROVER on custom data
python scripts/transfer_learning.py \
    --model grover \
    --data small_dataset.csv \
    --target activity \
    --task-type classification \
    --epochs 20
```

## 常见模式和最佳实践

### 模式 1：始终对分子使用支架分裂
```python
# GOOD: Prevents data leakage
splitter = dc.splits.ScaffoldSplitter()
train, test = splitter.train_test_split(dataset)

# BAD: Similar molecules in train and test
splitter = dc.splits.RandomSplitter()
train, test = splitter.train_test_split(dataset)
```

### 模式 2：标准化特征和目标
```python
transformers = [
    dc.trans.NormalizationTransformer(
        transform_y=True,  # Also normalize target values
        dataset=train
    )
]
for transformer in transformers:
    train = transformer.transform(train)
    test = transformer.transform(test)
```

### 模式 3：从简单开始，然后扩展
1.从随机森林+CircularFingerprint开始（快速基线）
2.如果RF效果良好，尝试XGBoost/LightGBM
3. 如果您有超过 5K 的样本，请转向深度学习（多任务回归器）
4. 如果你有超过 10K 的样本，请尝试 GNN
5. 对小型数据集或新颖的支架使用迁移学习

### 模式 4：处理不平衡数据
```python
# Option 1: Balancing transformer
transformer = dc.trans.BalancingTransformer(dataset=train)
train = transformer.transform(train)

# Option 2: Use balanced metrics
metric = dc.metrics.Metric(dc.metrics.balanced_accuracy_score)
```

### 模式 5：避免内存问题
```python
# Use DiskDataset for large datasets
dataset = dc.data.DiskDataset.from_numpy(X, y, w, ids)

# Use smaller batch sizes
model = dc.models.GCNModel(batch_size=32)  # Instead of 128
```

## 常见陷阱

### 问题 1：药物发现中的数据泄露
**问题**：使用随机分割允许训练/测试集中出现相似的分子。
**解决方案**：对于分子数据集始终使用 `ScaffoldSplitter`。

### 问题 2：GNN 与指纹相比表现不佳
**问题**：图神经网络的性能比简单指纹差。
**解决方案**：
- 确保数据集足够大（通常> 10K 样本）
- 增加训练次数（50-100）
- 尝试不同的架构（AttentiveFP、DMPNN 而不是 GCN）
- 使用预训练模型（GROVER）

### 问题 3：小数据集上的过度拟合
**问题**：模型记住训练数据。
**解决方案**：
- 使用更强的正则化（将 dropout 增加到 0.5）
- 使用更简单的模型（随机森林而不是深度学习）
- 应用迁移学习（ChemBERTa、GROVER）
- 收集更多数据

### 问题 4：导入错误
**问题**：模块未找到错误。
**解决方案**：确保 DeepChem 安装了所需的依赖项：
```bash
uv pip install deepchem
# For PyTorch models
uv pip install deepchem[torch]
# For all features
uv pip install deepchem[all]
```

## 参考文档

该技能包括全面的参考文档：

### `references/api_reference.md`
完整的 API 文档包括：
- 所有数据加载器及其用例
- 数据集类别以及何时使用每个类别
- 完整的特征器目录和选择指南
- 按类别组织的模型目录（50 多个模型）
- MoleculeNet 数据集描述
- 指标和评估函数
- 常见的代码模式

**何时参考**：当您需要特定 API 详细信息、参数名称或想要探索可用选项时，搜索此文件。

### `references/workflows.md`
八个详细的端到端工作流程：
1. SMILES 的分子性质预测
2. 使用 MoleculeNet 基准测试
3. 超参数优化
4. 使用预训练模型进行迁移学习
5. GAN 的分子生成
6. 材料性能预测
7. 蛋白质序列分析
8. 自定义模型集成

**何时参考**：使用这些工作流程作为实施完整解决方案的模板。

## 安装注意事项

基本安装：
```bash
uv pip install deepchem
```

对于 PyTorch 模型（GCN、GAT 等）：
```bash
uv pip install deepchem[torch]
```

对于所有功能：
```bash
uv pip install deepchem[all]
```

如果发生导入错误，用户可能需要特定的依赖项。查看 DeepChem 文档以获取详细的安装说明。

## 其他资源

- 官方文档：https://deepchem.readthedocs.io/
- GitHub 存储库：https://github.com/deepchem/deepchem
- 教程：https://deepchem.readthedocs.io/en/latest/get_started/tutorials.html
- 论文：“MoleculeNet：分子机器学习的基准”