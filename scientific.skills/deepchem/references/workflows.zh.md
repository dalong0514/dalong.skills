<!-- 此文件由机器翻译自 workflows.md -->

# DeepChem 工作流程

本文档提供常见 DeepChem 使用案例的详细工作流程。

## 工作流程 1：从 SMILES 预测分子特性

**目标**：根据 SMILES 字符串预测分子特性（例如溶解度、毒性、活性）。

### 分步过程

#### 1. 准备您的数据
数据应采用 CSV 格式，至少包含：
- 带有 SMILES 字符串的列
- 具有属性值（目标）的一列或多列

CSV 结构示例：
```csv
smiles,solubility,toxicity
CCO,-0.77,0
CC(=O)OC1=CC=CC=C1C(=O)O,-1.19,1
```

#### 2. 选择特征器
决策树：
- **小数据集 (<1K)**：使用 `CircularFingerprint` 或 `RDKitDescriptors`
- **中型数据集 (1K-100K)**：使用 `CircularFingerprint` 或 `MolGraphConvFeaturizer`
- **大型数据集 (>100K)**：使用基于图形的特征化器 (`MolGraphConvFeaturizer`、`DMPNNFeaturizer`)
- **迁移学习**：使用预训练的模型特征器 (`GroverFeaturizer`)

#### 3. 加载数据并对其进行特征化
<<<代码块_1>>>

#### 4. 分割数据
**关键**：使用 `ScaffoldSplitter` 进行药物发现以防止数据泄露。

<<<代码块_2>>>

#### 5. 转换数据（可选但推荐）
<<<代码块_3>>>

#### 6. 选择并训练模型
<<<代码块_4>>>

#### 7. 评估
<<<代码块_5>>>

#### 8. 做出预测
<<<代码块_6>>>

---

## 工作流程 2：使用 MoleculeNet 基准数据集

**目标**：在标准基准上快速训练和评估模型。

### 快速入门
```python
import deepchem as dc

# Load benchmark dataset
tasks, datasets, transformers = dc.molnet.load_tox21(
    featurizer='GraphConv',
    splitter='scaffold'
)
train, valid, test = datasets

# Train model
model = dc.models.GCNModel(
    n_tasks=len(tasks),
    mode='classification'
)
model.fit(train, nb_epoch=50)

# Evaluate
metric = dc.metrics.Metric(dc.metrics.roc_auc_score)
test_score = model.evaluate(test, [metric])
print(f"Test ROC-AUC: {test_score}")
```

### 可用的特征化器选项
调用 `load_*()` 函数时：
- `'ECFP'`：扩展连接指纹（圆形指纹）
- `'GraphConv'`：图卷积特征
- `'Weave'`：编织特征
- `'Raw'`：原始 SMILES 字符串
- `'smiles2img'`：2D 分子图像

### 可用的分离器选项
- `'scaffold'`：基于支架的拆分（推荐用于药物发现）
- `'random'`：随机分割
- `'stratified'`：分层分割（保留类别分布）
- `'butina'`：基于 Butina 聚类的分割

---

## 工作流程 3：超参数优化

**目标**：系统地找到最佳模型超参数。

### 使用 GridHyperparamOpt
```python
import deepchem as dc
import numpy as np

# Load data
tasks, datasets, transformers = dc.molnet.load_bbbp(
    featurizer='ECFP',
    splitter='scaffold'
)
train, valid, test = datasets

# Define parameter grid
params_dict = {
    'layer_sizes': [[1000], [1000, 500], [1000, 1000]],
    'dropouts': [0.0, 0.25, 0.5],
    'learning_rate': [0.001, 0.0001]
}

# Define model builder function
def model_builder(model_params, model_dir):
    return dc.models.MultitaskClassifier(
        n_tasks=len(tasks),
        n_features=1024,
        **model_params
    )

# Setup optimizer
metric = dc.metrics.Metric(dc.metrics.roc_auc_score)
optimizer = dc.hyper.GridHyperparamOpt(model_builder)

# Run optimization
best_model, best_params, all_results = optimizer.hyperparam_search(
    params_dict,
    train,
    valid,
    metric,
    transformers=transformers
)

print(f"Best parameters: {best_params}")
print(f"Best validation score: {all_results['best_validation_score']}")
```

---

## 工作流程 4：使用预训练模型进行迁移学习

**目标**：利用预训练模型来提高小数据集上的性能。

### 使用 ChemBERTa
```python
import deepchem as dc
from transformers import AutoTokenizer

# Load your data
loader = dc.data.CSVLoader(
    tasks=['activity'],
    feature_field='smiles',
    featurizer=dc.feat.DummyFeaturizer()  # ChemBERTa handles featurization
)
dataset = loader.create_dataset('data.csv')

# Split data
splitter = dc.splits.ScaffoldSplitter()
train, test = splitter.train_test_split(dataset)

# Load pretrained ChemBERTa
model = dc.models.HuggingFaceModel(
    model='seyonec/ChemBERTa-zinc-base-v1',
    task='regression',
    n_tasks=1
)

# Fine-tune
model.fit(train, nb_epoch=10)

# Evaluate
predictions = model.predict(test)
```

### 使用 GROVER
```python
# GROVER: pre-trained on molecular graphs
model = dc.models.GroverModel(
    task='classification',
    n_tasks=1,
    model_dir='./grover_model'
)

# Fine-tune on your data
model.fit(train_dataset, nb_epoch=20)
```

---

## 工作流程 5：使用 GAN 进行分子生成

**目标**：生成具有所需特性的新型分子。

### 基本 MolGAN
```python
import deepchem as dc

# Load training data (molecules for the generator to learn from)
tasks, datasets, _ = dc.molnet.load_qm9(
    featurizer='GraphConv',
    splitter='random'
)
train, _, _ = datasets

# Create and train MolGAN
gan = dc.models.BasicMolGANModel(
    learning_rate=0.001,
    vertices=9,  # max atoms in molecule
    edges=5,     # max bonds
    nodes=[128, 256, 512]
)

# Train
gan.fit_gan(
    train,
    nb_epoch=100,
    generator_steps=0.2,
    checkpoint_interval=10
)

# Generate new molecules
generated_molecules = gan.predict_gan_generator(1000)
```

### 条件生成
```python
# For property-targeted generation
from deepchem.models.optimizers import ExponentialDecay

gan = dc.models.BasicMolGANModel(
    learning_rate=ExponentialDecay(0.001, 0.9, 1000),
    conditional=True  # enable conditional generation
)

# Train with properties
gan.fit_gan(train, nb_epoch=100)

# Generate molecules with target properties
target_properties = np.array([[5.0, 300.0]])  # e.g., [logP, MW]
molecules = gan.predict_gan_generator(
    1000,
    conditional_inputs=target_properties
)
```

---

## 工作流程 6：材料特性预测

**目标**：预测晶体材料的特性。

### 使用 Crystal Graph 卷积网络
```python
import deepchem as dc

# Load materials data (structure files in CIF format)
loader = dc.data.CIFLoader()
dataset = loader.create_dataset('materials.csv')

# Split data
splitter = dc.splits.RandomSplitter()
train, test = splitter.train_test_split(dataset)

# Create CGCNN model
model = dc.models.CGCNNModel(
    n_tasks=1,
    mode='regression',
    batch_size=32,
    learning_rate=0.001
)

# Train
model.fit(train, nb_epoch=100)

# Evaluate
metric = dc.metrics.Metric(dc.metrics.mae_score)
test_score = model.evaluate(test, [metric])
```

---

## 工作流程 7：蛋白质序列分析

**目标**：从序列预测蛋白质特性。

### 使用 ProtBERT
```python
import deepchem as dc

# Load protein sequence data
loader = dc.data.FASTALoader()
dataset = loader.create_dataset('proteins.fasta')

# Use ProtBERT
model = dc.models.HuggingFaceModel(
    model='Rostlab/prot_bert',
    task='classification',
    n_tasks=1
)

# Split and train
splitter = dc.splits.RandomSplitter()
train, test = splitter.train_test_split(dataset)
model.fit(train, nb_epoch=5)

# Predict
predictions = model.predict(test)
```

---

## 工作流程 8：自定义模型集成

**目标**：将您自己的 PyTorch/scikit-learn 模型与 DeepChem 结合使用。

### 包装 Scikit-Learn 模型
```python
from sklearn.ensemble import RandomForestRegressor
import deepchem as dc

# Create scikit-learn model
sklearn_model = RandomForestRegressor(
    n_estimators=100,
    max_depth=10,
    random_state=42
)

# Wrap in DeepChem
model = dc.models.SklearnModel(model=sklearn_model)

# Use with DeepChem datasets
model.fit(train)
predictions = model.predict(test)

# Evaluate
metric = dc.metrics.Metric(dc.metrics.r2_score)
score = model.evaluate(test, [metric])
```

### 创建自定义 PyTorch 模型
```python
import torch
import torch.nn as nn
import deepchem as dc

class CustomNetwork(nn.Module):
    def __init__(self, n_features, n_tasks):
        super().__init__()
        self.fc1 = nn.Linear(n_features, 512)
        self.fc2 = nn.Linear(512, 256)
        self.fc3 = nn.Linear(256, n_tasks)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(0.2)

    def forward(self, x):
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.relu(self.fc2(x))
        x = self.dropout(x)
        return self.fc3(x)

# Wrap in DeepChem TorchModel
model = dc.models.TorchModel(
    model=CustomNetwork(n_features=2048, n_tasks=1),
    loss=nn.MSELoss(),
    output_types=['prediction']
)

# Train
model.fit(train, nb_epoch=50)
```

---

## 常见陷阱和解决方案

### 问题 1：药物发现中的数据泄露
**问题**：使用随机分割允许训练和测试集中出现相似的分子。
**解决方案**：对于分子数据集始终使用 `ScaffoldSplitter`。

### 问题 2：分类不平衡
**问题**：少数族裔班级表现不佳。
**解决方案**：使用 `BalancingTransformer` 或加权指标。
```python
transformer = dc.trans.BalancingTransformer(dataset=train)
train = transformer.transform(train)
```

### 问题 3：大型数据集的内存问题
**问题**：数据集不适合内存。
**解决方案**：使用 `DiskDataset` 而不是 `NumpyDataset`。
```python
dataset = dc.data.DiskDataset.from_numpy(X, y, w, ids)
```

### 问题 4：小数据集上的过度拟合
**问题**：模型记住训练数据。
**解决方案**：
1.使用更强的正则化（增加dropout）
2.使用更简单的模型（随机森林、岭）
3.应用迁移学习（预训练模型）
4.收集更多数据

### 问题 5：图神经网络性能不佳
**问题**：GNN 的性能比指纹差。
**解决方案**：
1. 检查数据集是否足够大（GNN 通常需要 >10K 样本）
2.增加训练次数
3.尝试不同的GNN架构（AttentiveFP、DMPNN）
4.使用预训练模型（GROVER）