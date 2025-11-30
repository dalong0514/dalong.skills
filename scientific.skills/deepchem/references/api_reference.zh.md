<!-- 此文件由机器翻译自 api_reference.md -->

# DeepChem API 参考

本文档提供了按功能组织的 DeepChem 核心 API 的综合参考。

## 数据处理

### 数据加载器

#### 文件格式加载器
- **CSVLoader**：通过可自定义的功能处理从 CSV 文件加载表格数据
- **UserCSVLoader**：用户定义的 CSV 加载，具有灵活的列规格
- **SDFLoader**：处理分子结构文件（SDF格式）
- **JsonLoader**：导入 JSON 结构的数据集
- **ImageLoader**：加载计算机视觉任务的图像数据

#### 生物数据加载器
- **FASTALoader**：处理 FASTA 格式的蛋白质/DNA 序列
- **FASTQLoader**：处理带有质量分数的 FASTQ 测序数据
- **SAMLoader/BAMLoader/CRAMLoader**：支持序列比对格式

#### 专用装载机
- **DFTYamlLoader**：处理密度泛函理论计算数据
- **InMemoryLoader**：直接从Python对象加载数据

### 数据集类

- **NumpyDataset**：包装 NumPy 数组以进行内存数据操作
- **DiskDataset**：管理存储在磁盘上的较大数据集，减少内存开销
- **ImageDataset**：用于基于图像的 ML 任务的专用容器

### 数据分割器

#### 通用分配器
- **RandomSplitter**：随机数据集分区
- **IndexSplitter**：按指定索引分割
- **SpecifiedSplitter**：使用预定义的分割
- **RandomStratifiedSplitter**：分层随机分割
- **SingletaskStratifiedSplitter**：单个任务的分层拆分
- **TaskSplitter**：针对多任务场景进行拆分

#### 分子特异性分离器
- **ScaffoldSplitter**：按结构支架划分分子（防止数据泄漏）
- **ButinaSplitter**：基于聚类的分子分裂
- **FingerprintSplitter**：根据分子指纹相似性进行分割
- **MaxMinSplitter**：最大化训练/测试集之间的多样性
- **MolecularWeightSplitter**：按分子量属性拆分

**最佳实践**：对于药物发现任务，使用 ScaffoldSplitter 来防止相似分子结构的过度拟合。

### 变形金刚

#### 标准化
- **NormalizationTransformer**：标准标准化（mean=0，std=1）
- **MinMaxTransformer**：将特征缩放到 [0,1] 范围
- **LogTransformer**：应用日志转换
- **PowerTransformer**：Box-Cox 和 Yeo-Johnson 变换
- **CDFTransformer**：累积分布函数标准化

#### 特定任务
- **BalancingTransformer**：解决类不平衡问题
- **FeaturizationTransformer**：应用动态特征工程
- **CoulombFitTransformer**：量子化学特定
- **DAGTransformer**：有向非循环图转换
- **RxnSplitTransformer**：化学反应预处理

## 分子特征分析

### 基于图形的特征器
将它们与图神经网络（GCN、MPNN 等）一起使用：

- **ConvMolFeaturizer**：图卷积网络的图表示
- **WeaveFeaturizer**：“Weave”图嵌入
- **MolGraphConvFeaturizer**：图卷积就绪表示
- **EquivariantGraphFeaturizer**：保持几何不变性
- **DMPNNFeaturizer**：定向消息传递神经网络输入
- **GroverFeaturizer**：预训练的分子嵌入

### 基于指纹的特征器
将它们与传统 ML（随机森林、SVM、XGBoost）结合使用：

- **MACCSKeysFingerprint**：167 位结构密钥
- **CircularFingerprint**：扩展连接指纹（摩根指纹）
  - 参数：`radius`（默认 2）、`size`（默认 2048）、`useChirality`（默认 False）
- **PubChemFingerprint**：881 位结构描述符
- **Mol2VecFingerprint**：学习的分子向量表示

### 描述符特征器
直接计算分子性质：

- **RDKitDescriptors**：约 200 个分子描述符（MW、LogP、H 供体、H 受体、TPSA 等）
- **MordredDescriptors**：综合结构和物理化学描述符
- **库仑矩阵**：3D 结构的原子间距离矩阵

### 基于序列的特征器
对于循环网络和变压器：

- **SmilesToSeq**：将 SMILES 字符串转换为序列
- **SmilesToImage**：从 SMILES 生成 2D 图像表示
- **RawFeaturizer**：不变地传递原始分子数据

### 选型指南

|使用案例|推荐特征器 |型号类型 |
|----------|------------------------|------------|
|图神经网络 | ConvMolFeaturizer、MolGraphConvFeaturizer | GCN、MPNN、GAT |
|传统机器学习 | CircularFingerprint、RDKitDescriptors |随机森林、XGBoost、SVM |
|深度学习（非图）|圆形指纹、Mol2VecFingerprint |密集网络，CNN |
|序列模型|微笑到序列 | LSTM、GRU、变压器 |
| 3D 分子结构 |库仑矩阵|专业 3D 模型 |
|快速基线| RDKit 描述符 |线性、山脊、套索 |

## 型号

### Scikit-Learn 集成
- **SklearnModel**：任何 scikit-learn 算法的包装器
  - 用法：`SklearnModel(model=RandomForestRegressor())`

### 梯度提升
- **GBDTModel**：梯度提升决策树（XGBoost、LightGBM）

### PyTorch 模型

#### 分子性质预测
- **MultitaskRegressor**：具有共享表示的多任务回归
- **MultitaskClassifier**：多任务分类
- **MultitaskFitTransformRegressor**：具有学习转换的回归
- **GCNModel**：图卷积网络
- **GATModel**：图注意力网络
- **AttentiveFPModel**：注意力指纹网络
- **DMPNNModel**：定向消息传递神经网络
- **GroverModel**：GROVER 预训练变压器
- **MATModel**：分子注意力转换器

#### 材料科学
- **CGCNN模型**：水晶图卷积网络
- **MEGNetModel**：材料图网络
- **LCNNModel**：用于材料的 Lattice CNN

#### 生成模型
- **GANModel**：生成对抗网络
- **WGAN模型**：Wasserstein GAN
- **BasicMolGANModel**：分子 GAN
- **LSTMGenerator**：基于 LSTM 的分子生成
- **SeqToSeqModel**：序列到序列模型

#### 物理模型
- **PINNModel**：基于物理的神经网络
- **HNNModel**：哈密顿神经网络
- **LNN**：拉格朗日神经网络
- **FNOModel**：傅立叶神经算子

#### 计算机视觉
- **CNN**：卷积神经网络
- **UNetModel**：用于分割的 U-Net 架构
- **InceptionV3Model**：预训练的 Inception v3
- **MobileNetV2Model**：轻量级移动网络

### 拥抱脸部模特

- **HuggingFaceModel**：高频变压器的通用包装
- **Chemberta**：用于分子特性预测的化学 BERT
- **MoLFormer**：分子变压器架构
- **ProtBERT**：蛋白质序列 BERT
- **DeepAbLLM**：抗体大语言模型

### 选型指南

|任务|推荐型号 |特征化器 |
|------|--------------------|------------|
|小数据集（<1000 个样本）| SklearnModel（随机森林）|圆形指纹|
|中等数据集 (1K-100K) | GBDT模型或多任务回归器| CircularFingerprint 或 ConvMolFeaturizer |
|大型数据集 (>100K) | GCN 模型、AttentiveFP 模型或 DMPNN | MolGraphConvFeaturizer | MolGraphConvFeaturizer |
|迁移学习 | GroverModel，Chemberta，MoL前任 |特定型号|
|材料特性| CGCNN模型、MEGNet模型 |基于结构|
|分子生成| BasicMolGAN 模型、LSTM 生成器 |微笑到序列 |
|蛋白质序列|普特伯特 |基于序列|

## MoleculeNet 数据集

通过 `dc.molnet.load_*()` 函数快速访问 30 多个基准数据集。

### 分类数据集
- **load_bace()**：BACE-1 抑制剂（二元分类）
- **load_bbbp()**：血脑屏障穿透
- **load_clintox()**：临床毒性
- **load_hiv()**：HIV 抑制活性
- **load_muv()**：PubChem BioAssay（具有挑战性，稀疏）
- **load_pcba()**：PubChem 筛选数据
- **load_sider()**：药物不良反应（多标签）
- **load_tox21()**：12 种毒性测定（多任务）
- **load_toxcast()**：EPA ToxCast 筛查

### 回归数据集
- **load_delaney()**：水溶性 (ESOL)
- **load_freesolv()**：溶剂化自由能
- **load_lipo()**：亲脂性（辛醇-水分配）
- **load_qm7/qm8/qm9()**：量子力学性能
- **load_hopv()**：有机光伏特性

### 蛋白质-配体结合
- **load_pdbbind()**：绑定亲和力数据

### 材料科学
- **load_perovskite()**：钙钛矿稳定性
- **load_mp_formation_energy()**：材料项目形成能量
- **load_mp_metalicity()**：金属与非金属分类
- **load_bandgap()**：电子带隙预测

### 化学反应
- **load_uspto()**：USPTO 反应数据集

### 使用模式
```python
tasks, datasets, transformers = dc.molnet.load_bbbp(
    featurizer='GraphConv',  # or 'ECFP', 'GraphConv', 'Weave', etc.
    splitter='scaffold',      # or 'random', 'stratified', etc.
    reload=False              # set True to skip caching
)
train, valid, test = datasets
```

## 指标

`dc.metrics` 中提供的常见评估指标：

### 分类指标
- **roc_auc_score**：ROC 曲线下的面积（二元/多类）
- **prc_auc_score**：精确率-召回率曲线下的面积
- **accuracy_score**：分类准确度
- **balanced_accuracy_score**：不平衡数据集的平衡准确性
- **recall_score**：灵敏度/召回率
- **精度_分数**：精度
- **f1_score**：F1 分数

### 回归指标
- **平均绝对误差**：MAE
- **均方误差**：MSE
- **均方根误差**：RMSE
- **r2_score**：R² 决定系数
- **pearson_r2_score**：皮尔逊相关性
- **spearman_correlation**：斯皮尔曼等级相关

### 多任务指标
大多数指标通过对任务进行平均来支持多任务评估。

## 训练模式

标准 DeepChem 工作流程：

<<<代码块_1>>>

## 常见模式

### 模式 1：使用 MoleculeNet 进行快速基线
<<<代码块_2>>>

### 模式 2：使用图网络的自定义数据
<<<代码块_3>>>

### 模式 3：使用预训练模型进行迁移学习
<<<代码块_4>>>