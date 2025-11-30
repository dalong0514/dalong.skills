<!-- 此文件由机器翻译自 api_reference.md -->

# Molfeat API 参考

## 核心模块

Molfeat 分为几个关键模块，提供分子特征化的不同方面：

- **`molfeat.store`** - 管理模型加载、列出和注册
- **`molfeat.calc`** - 提供单分子特征化计算器
- **`molfeat.trans`** - 提供用于批处理的 scikit-learn 兼容变压器
- **`molfeat.utils`** - 用于数据处理的实用函数
- **`molfeat.viz`>** - 分子特征的可视化工具

---

## molfeat.calc - 计算器

计算器是可调用对象，可将单个分子转换为特征向量。它们接受 RDKit `Chem.Mol` 对象或 SMILES 字符串作为输入。

### SerializedCalculator（基类）

所有计算器的基本抽象类。子类化时，必须实现：
- `__call__()` - 特征化所需的方法
- `__len__()` - 可选，返回输出长度
- `columns` - 可选属性，返回功能名称
- `batch_compute()` - 可选，用于高效的批处理

**状态管理方法：**
- `to_state_json()` - 将计算器状态保存为 JSON
- `to_state_yaml()` - 将计算器状态保存为 YAML
- `from_state_dict()` - 从状态字典加载计算器
- `to_state_dict()` - 将计算器状态导出为字典

### FP计算器

计算分子指纹。支持15+指纹方式。

**支持的指纹类型：**

**结构指纹：**
- `ecfp` - 扩展连接指纹（圆形）
- `fcfp` - 功能级指纹
- `rdkit` - RDKit 拓扑指纹
- `maccs` - MACCS 密钥（166 位结构密钥）
- `avalon` - Avalon 指纹
- `pattern` - 模式指纹
- `layered` - 分层指纹

**基于原子的指纹：**
- `atompair` - 原子对指纹
- `atompair-count` - 计算原子对
- `topological` - 拓扑扭转指纹
- `topological-count` - 计算拓扑扭转

**专业指纹：**
- `map4` - MinHashed 原子对指纹最多 4 个键
- `secfp` - SMILES 扩展连接指纹
- `erg` - 扩展简化图
- `estate` - 电拓扑状态索引

**参数：**
- `method` (str) - 指纹类型名称
- `radius` (int) - 圆形指纹的半径（默认值：3）
- `fpSize` (int) - 指纹大小（默认值：2048）
- `includeChirality` (bool) - 包括手性信息
- `counting` (bool) - 使用计数向量而不是二进制

**用途：**
```python
from molfeat.calc import FPCalculator

# Create fingerprint calculator
calc = FPCalculator("ecfp", radius=3, fpSize=2048)

# Compute fingerprint for single molecule
fp = calc("CCO")  # Returns numpy array

# Get fingerprint length
length = len(calc)  # 2048

# Get feature names
names = calc.columns
```

**常见指纹尺寸：**
- MACCS：167 个维度
- ECFP（默认）：2048 维
- MAP4（默认）：1024 维

### 描述符计算器

**RDKitDescriptors2D**
使用 RDKit 计算 2D 分子描述符。

<<<代码块_1>>>

**RDKitDescriptors3D**
计算 3D 分子描述符（需要构象异构体生成）。

**莫德雷德描述符**
使用 Mordred 计算超过 1800 个分子描述符。

<<<代码块_2>>>

### 药效团计算器

**药效团二维**
RDKit 的 2D 药效团指纹生成。

**药效基团3D**
来自多个构象异构体的共有药效团指纹。

**CATS计算器**
计算化学高级模板搜索 (CATS) 描述符 - 药效团点对分布。

**参数：**
- `mode` - “2D”或“3D”距离计算
- `dist_bins` - 对分布的距离箱
- `scale` - 缩放模式：“raw”、“num”或“count”

<<<代码块_3>>>

### 形状描述符

**USR描述符**
超快形状识别描述符（多个变体）。

**ElectroShape描述符**
静电形状描述符结合了形状、手性和静电学。

### 基于图形的计算器

**脚手架钥匙计算器**
计算 40 多种基于支架的分子特性。

**原子计算器**
图神经网络的原子级特征化。

**债券计算器**
图神经网络的键级特征化。

### 实用函数

**获取计算器()**
通过名称实例化计算器的工厂函数。

<<<代码块_4>>>

对于不支持的特征器，提高 `ValueError`。

---

## molfeat.trans - 变形金刚
Transformer 将计算器包装到完整的特征化管道中以进行批处理。

### 分子转换器

用于批量分子特征化的 Scikit-learn 兼容转换器。

**关键参数：**
- `featurizer` - 使用的计算器或特征器
- `n_jobs` (int) - 并行作业数量（所有核心为-1）
- `dtype` - 输出数据类型（numpy float32/64，torch 张量）
- `verbose` (bool) - 启用详细日志记录
- `ignore_errors` (bool) - 失败时继续（对于失败的分子返回 None）

**基本方法：**
- `transform(mols)` - 处理批次并返回表示
- `_transform(mol)` - 处理单个分子特征化
- `__call__(mols)` - Transform() 的便捷包装
- `preprocess(mol)` - 准备输入分子（不自动应用）
- `to_state_yaml_file(path)` - 保存变压器配置
- `from_state_yaml_file(path)` - 加载变压器配置

**用途：**
<<<代码块_5>>>

**性能：** 对 642 个分子进行的测试显示，使用 4 个并行作业与单线程处理相比，速度提高了 3.4 倍。

### 壮举Concat

将多个特征化器连接成统一的表示。

<<<代码块_6>>>

### PretrainedMolTransformer

用于预训练深度学习模型的 `MoleculeTransformer` 的子类。

**独特的特点：**
- `_embed()` - 神经网络的批量推理
- `_convert()` - 将 SMILES/分子转换为模型兼容的格式
  - 用于语言模型的 SELFIES 字符串
  - 图神经网络的 DGL 图
- 集成缓存系统，实现高效存储

**用途：**
```python
from molfeat.trans.pretrained import PretrainedMolTransformer

# Load pretrained model
transformer = PretrainedMolTransformer("ChemBERTa-77M-MLM", n_jobs=-1)

# Generate embeddings
embeddings = transformer(smiles)
```

### 预计算MolTransformer

用于缓存/预计算功能的转换器。

---

## molfeat.store - 模型商店

管理特征器发现、加载和注册。

### 模型商店

用于访问可用特征器的中央枢纽。

**关键方法：**
- `available_models` - 列出所有可用特征的属性
- `search(name=None, **kwargs)` - 搜索特定特征器
- `load(name, **kwargs)` - 按名称加载特征器
- `register(name, card)` - 注册自定义特征器

**用途：**
```python
from molfeat.store.modelstore import ModelStore

# Initialize store
store = ModelStore()

# List all available models
all_models = store.available_models
print(f"Found {len(all_models)} featurizers")

# Search for specific model
results = store.search(name="ChemBERTa-77M-MLM")
if results:
    model_card = results[0]

    # View usage information
    model_card.usage()

    # Load the model
    transformer = model_card.load()

# Direct loading
transformer = store.load("ChemBERTa-77M-MLM")
```

**模型卡属性：**
- `name` - 型号标识符
- `description` - 型号描述
- `version` - 模型版本
- `authors` - 模型作者
- `tags` - 分类标签
- `usage()` - 显示使用示例
- `load(**kwargs)` - 加载模型

---

## 常见模式

### 错误处理

```python
# Enable error tolerance
featurizer = MoleculeTransformer(
    calc,
    n_jobs=-1,
    verbose=True,
    ignore_errors=True
)

# Failed molecules return None
features = featurizer(smiles_with_errors)
```

### 数据类型控制

```python
# NumPy float32 (default)
features = transformer(smiles, enforce_dtype=True)

# PyTorch tensors
import torch
transformer = MoleculeTransformer(calc, dtype=torch.float32)
features = transformer(smiles)
```

### 持久性和再现性

```python
# Save transformer state
transformer.to_state_yaml_file("config.yml")
transformer.to_state_json_file("config.json")

# Load from saved state
transformer = MoleculeTransformer.from_state_yaml_file("config.yml")
transformer = MoleculeTransformer.from_state_json_file("config.json")
```

### 预处理

```python
# Manual preprocessing
mol = transformer.preprocess("CCO")

# Transform with preprocessing
features = transformer.transform(smiles_list)
```

---

## 集成示例

### Scikit-learn 管道

```python
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from molfeat.trans import MoleculeTransformer
from molfeat.calc import FPCalculator

# Create pipeline
pipeline = Pipeline([
    ('featurizer', MoleculeTransformer(FPCalculator("ecfp"))),
    ('classifier', RandomForestClassifier())
])

# Fit and predict
pipeline.fit(smiles_train, y_train)
predictions = pipeline.predict(smiles_test)
```

### PyTorch 集成

```python
import torch
from torch.utils.data import Dataset, DataLoader
from molfeat.trans import MoleculeTransformer

class MoleculeDataset(Dataset):
    def __init__(self, smiles, labels, transformer):
        self.smiles = smiles
        self.labels = labels
        self.transformer = transformer

    def __len__(self):
        return len(self.smiles)

    def __getitem__(self, idx):
        features = self.transformer(self.smiles[idx])
        return torch.tensor(features), torch.tensor(self.labels[idx])

# Create dataset and dataloader
transformer = MoleculeTransformer(FPCalculator("ecfp"))
dataset = MoleculeDataset(smiles, labels, transformer)
loader = DataLoader(dataset, batch_size=32)
```

---

## 性能提示

1. **并行化**：使用`n_jobs=-1`来利用所有CPU核心
2. **批处理**：一次处理多个分子而不是循环
3. **缓存**：利用预训练模型的内置缓存
4. **数据类型**：在精度允许的情况下使用 float32 而不是 float64
5. **错误处理**：为具有潜在无效分子的大型数据集设置`ignore_errors=True`