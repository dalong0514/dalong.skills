<!-- 此文件由机器翻译自 available_featurizers.md -->

# Molfeat 中可用的特征器

本文档提供了 molfeat 中所有可用特征器的综合目录，按类别组织。

## 基于 Transformer 的语言模型

使用 SMILES/SELFIES 表示形式进行分子嵌入的预训练 Transformer 模型。

### RoBERTa 风格的模型
- **Roberta-Zinc480M-102M** - RoBERTa 掩码语言模型在 ZINC 数据库中的约 480M SMILES 字符串上进行训练
- **ChemBERTa-77M-MLM** - 基于 RoBERTa 的屏蔽语言模型，在 77M PubChem 化合物上进行训练
- **ChemBERTa-77M-MTR** - 在 PubChem 化合物上训练的多任务回归版本

### GPT 式自回归模型
- **GPT2-Zinc480M-87M** - 在 ZINC 的约 480M SMILES 上训练的 GPT-2 自回归语言模型
- **ChemGPT-1.2B** - 在 PubChem10M 上预训练的大型变压器（1.2B 参数）
- **ChemGPT-19M** - 在 PubChem10M 上预训练的中型变压器（19M 参数）
- **ChemGPT-4.7M** - 在 PubChem10M 上预训练的小型变压器（4.7M 参数）

### 专用变压器模型
- **MolT5** - 用于分子字幕和基于文本的生成的自监督框架

## 图神经网络 (GNN)

在分子图结构上运行的预训练图神经网络模型。

### GIN（图同构网络）变体
所有这些都针对具有不同目标的 ChEMBL 分子进行了预训练：
- **gin-supervised-masking** - 使用节点屏蔽目标进行监督
- **gin-supervised-infomax** - 通过图级互信息最大化进行监督
- **gin-supervised-edgepred** - 使用边缘预测目标进行监督
- **gin-supervised-contextpred** - 使用上下文预测目标进行监督

### 其他基于图的模型
- **JTVAE_zinc_no_kl** - 用于分子生成的连接树 VAE（在 ZINC 上训练）
- **Graphormer-pcqm4mv2** - 在 PCQM4Mv2 量子化学数据集上预训练的图形转换器，用于 HOMO-LUMO 间隙预测

## 分子描述符

物理化学性质和分子特征计算器。

### 2D 描述符
- **desc2D** / **rdkit2D** - 200+ RDKit 2D 分子描述符，包括：
  - 分子量、logP、TPSA
  - H键供体/受体
  - 可旋转债券
  - 环数和芳香度
  - 分子复杂性指标

### 3D 描述符
- **desc3D** / **rdkit3D** - RDKit 3D 分子描述符（需要构象异构体生成）
  - 惯性矩
  - PMI（主惯性矩）比率
  - 非球面、偏心率
  - 回转半径

### 综合描述符集
- **莫德雷德** - 超过 1800 个分子描述符，涵盖：
  - 宪法描述符
  - 拓扑索引
  - 连接性指数
  - 信息内容
  - 2D/3D 自相关
  - 突发奇想描述符
  - 逍遥游描述符
  - 还有更多

### 电拓扑描述符
- **estate** - 电拓扑状态（E-State）索引编码：
  - 原子环境信息
  - 电子和拓扑特性
  - 杂原子贡献

## 分子指纹

表示分子子结构的二进制或基于计数的固定长度向量。

### 圆形指纹（ECFP 式）
- **ecfp** / **ecfp:2** / **ecfp:4** / **ecfp:6** - 扩展连接指纹
  - 半径变体（2、4、6 对应于直径）
  - 默认：半径=3，2048 位
  - 最流行的相似性搜索
- **ecfp-count** - ECFP 的计数版本（非二进制）
- **fcfp** / **fcfp-count** - 功能级圆形指纹
  - 与 ECFP 类似，但使用功能组
  - 更好地基于药效基团的相似性

### 基于路径的指纹
- **rdkit** - 基于线性路径的RDKit拓扑指纹
- **模式** - 模式指纹（类似于 MACCS 但自动化）
- **分层** - 具有多个子结构层的分层指纹

### 基于密钥的指纹
- **maccs** - MACCS 密钥（166 位结构密钥）
  - 修复了一组预定义的子结构
  - 适合脚手架跳跃
  - 快速计算
- **阿瓦隆** - 阿瓦隆指纹
  - 与 MACCS 类似，但功能更多
  - 优化相似性搜索

### 原子对指纹
- **atompair** - 原子对指纹
  - 编码原子对和它们之间的距离
  - 适合 3D 相似度
- **atompair-count** - 原子对的计数版本

### 拓扑扭转指纹
- **拓扑** - 拓扑扭转指纹
  - 编码 4 个相连原子的序列
  - 捕获本地拓扑
- **topological-count** - 拓扑扭转的计数版本
### MinHashed 指纹
- **map4** - MinHashed 原子对指纹最多 4 个键
  - 结合了原子对和 ECFP 概念
  - 默认：1024 维
  - 快速高效地处理大型数据集
- **secfp** - SMILES 扩展连接指纹
  - 直接对 SMILES 字符串进行操作
  - 捕获子结构和原子对信息

### 扩展简化图
- **erg** - 扩展简化图
  - 使用药效点代替原子
  - 降低图形复杂性，同时保留关键特征

## 药效团描述符

基于药理学相关官能团及其空间关系的特征。

### CATS（化学高级模板搜索）
- **cats2D** - 2D CATS 描述符
  - 药效团点对分布
  - 基于最短路径的距离
  - 默认21个描述符
- **cats3D** - 3D CATS 描述符
  - 基于欧几里德距离
  - 需要构象异构体生成
- **cats2D_pharm** / **cats3D_pharm** - 药效团变体

### 戈比药效团
- **gobbi2D** - 2D 药效团指纹
  - 8 种药效基团特征类型：
    - 疏水性
    - 芳香
    - H键受体
    - H键供体
    - 可正电离
    - 负电离
    - 集总疏水物
  - 适合虚拟筛选

### Pmapper 药效团
- **pmapper2D** - 2D 药效团特征
- **pmapper3D** - 3D 药效团签名
  - 高维药效团描述符
  - 对于 QSAR 和相似性搜索有用

## 形状描述符

捕获 3D 分子形状和静电特性的描述符。

### USR（超快形状识别）
- **usr** - 基本 USR 描述符
  - 12维编码形状分布
  - 极快的计算
- **usrcat** - 具有药效团限制的 USR
  - 60 个维度（每个特征类型 12 个）
  - 结合形状和药效团信息

### 静电形状
- ** electroshape ** - ElectroShape 描述符
  - 结合了分子形状、手性和静电学
  - 对于蛋白质-配体对接预测有用

## 基于支架的描述符

基于分子支架和核心结构的描述符。

### 脚手架钥匙
- **脚手架钥匙** - 脚手架钥匙计算器
  - 40 多个基于脚手架的属性
  - 生物电子等排支架表示
  - 捕捉核心结构特征

## GNN 输入的图特征器

用于构建图神经网络的图表示的原子和键级特征。

### 原子级特征
- **atom-onehot** - One-hot 编码原子特征
- **atom-default** - 默认原子特征包括：
  - 原子序数
  - 学位、正式费用
  - 杂交
  - 芳香度
  - 氢原子数

### 债券级特征
- **bond-onehot** - One-hot 编码键功能
- **债券默认** - 默认债券特征包括：
  - 键类型（单键、双键、三键、芳香键）
  - 共轭
  - 戒指会员资格
  - 立体化学

## 集成预训练模型集合

Molfeat 集成了各种来源的模型：

### HuggingFace 模型
通过 HuggingFace hub 访问变压器模型：
- ChemBERTa 变体
- ChemGPT 变体
- 摩尔T5
- 自定义上传的模型

### DGL-LifeSci 模型
来自 DGL-Life 的预训练 GNN 模型：
- 具有不同预训练任务的 GIN 变体
- 细心的FP模型
- MPNN模型

### FCD（Fréchet ChemNet 距离）
- **fcd** - 用于分子生成评估的预训练 CNN

### Graphomer 模型
- 来自微软研究院的图形转换器
- 在量子化学数据集上进行预训练

## 使用说明

### 选择特征器

**对于传统机器学习（随机森林、SVM 等）：**
- 从 **ecfp** 或 **maccs** 指纹开始
- 尝试 **desc2D** 可解释模型
- 使用 **FeatConcat** 组合多个指纹

**对于深度学习：**
- 使用 **ChemBERTa** 或 **ChemGPT** 进行变压器嵌入
- 使用 **gin-supervised-*** 进行图神经网络嵌入
- 考虑使用 **Graphormer** 进行量子属性预测

**对于相似性搜索：**
- **ecfp** - 通用，最受欢迎
- **maccs** - 快速，适合支架跳跃
- **map4** - 高效进行大规模搜索
- **usr** / **usrcat** - 3D 形状相似度

**对于基于药效团的方法：**
- **fcfp** - 基于功能组
- **cats2D/3D** - 药效团对分布
- **gobbi2D** - 显式药效团特征

**为了可解释性：**
- **desc2D** / **mordred** - 命名描述符
- **maccs** - 可解释的子结构键
- **脚手架键** - 基于脚手架的功能

### 模型依赖关系

一些特征器需要可选的依赖项：

- **DGL 模型**（gin-*、jtvae）：`pip install "molfeat[dgl]"`
- **绘图者**：`pip install "molfeat[graphormer]"`
- **变形金刚**（ChemBERTa、ChemGPT、MolT5）：`pip install "molfeat[transformer]"`
- **FCD**：`pip install "molfeat[fcd]"`
- **MAP4**：`pip install "molfeat[map4]"`
- **所有依赖项**：`pip install "molfeat[all]"`

### 访问所有可用模型

```python
from molfeat.store.modelstore import ModelStore

store = ModelStore()
all_models = store.available_models

# Print all available featurizers
for model in all_models:
    print(f"{model.name}: {model.description}")

# Search for specific types
transformers = [m for m in all_models if "transformer" in m.tags]
gnn_models = [m for m in all_models if "gnn" in m.tags]
fingerprints = [m for m in all_models if "fingerprint" in m.tags]
```

## 性能特点

### 计算速度（相对）
**最快：**
- 麦克斯
- ECFP
- rdkit指纹
- 用户

**中：**
- 描述2D
- 猫2D
- 大多数指纹

**较慢：**
- 莫德雷德（1800+ 描述符）
- desc3D（需要生成构象异构体）
- 一般 3D 描述符

**最慢（第一次运行）：**
- 预训练模型（ChemBERTa、ChemGPT、GIN）
- 注意：后续运行受益于缓存

### 维度

**低（< 200 暗淡）：**
- mac (167)
- 用户 (12)
- 乌尔卡特 (60)

**中（200-2000 暗淡）：**
- desc2D (~200)
- ecfp（默认2048，可配置）
-map4（默认1024）

**高（> 2000 暗淡）：**
- 莫德雷德 (1800+)
- 串联指纹
- 一些变压器嵌入

**变量：**
- 变压器型号（通常为 768-1024）
- GNN 模型（取决于架构）