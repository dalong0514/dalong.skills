<!-- 此文件由机器翻译自 models.md -->

# PyHealth 模型

## 概述

PyHealth 为医疗保健预测任务提供了 33 多个模型，从简单的基线到最先进的深度学习架构。模型分为通用架构和医疗保健特定模型。

## 模型基类

所有模型均继承自 `BaseModel` 并具有标准 PyTorch 功能：

**关键属性：**
- `dataset`：关联的 SampleDataset
- `feature_keys`：要使用的输入功能（例如，[“诊断”、“药物”]）
- `mode`：任务类型（“二进制”、“多类”、“多标签”、“回归”）
- `embedding_dim`：特征嵌入维度
- `device`：计算设备（CPU/GPU）

**关键方法：**
- `forward()`：模型前向传递
- `train_step()`：单次训练迭代
- `eval_step()`：单次评估迭代
- `save()`：保存模型检查点
- `load()`：加载模型检查点

## 通用型号

### 基线模型

**逻辑回归** (`LogisticRegression`)
- 具有均值池化的线性分类器
- 简单的比较基线
- 快速训练和推理
- 有利于可解释性

**用途：**
```python
from pyhealth.models import LogisticRegression

model = LogisticRegression(
    dataset=sample_dataset,
    feature_keys=["diagnoses", "medications"],
    mode="binary"
)
```

**多层感知器** (`MLP`)
- 前馈神经网络
- 可配置的隐藏层
- 支持平均值/总和/最大池化
- 结构化数据的良好基线

**参数：**
- `hidden_dim`：隐藏层大小
- `num_layers`：隐藏层数
- `dropout`：辍学率
- `pooling`：聚合方法（“平均值”、“总和”、“最大值”）

**用途：**
<<<代码块_1>>>

### 卷积神经网络

**CNN** (`CNN`)
- 用于模式检测的卷积层
- 对顺序数据和空间数据有效
- 捕捉局部时间模式
- 参数高效

**架构：**
- 多个一维卷积层
- 最大池化以减少维度
- 全连接的输出层

**参数：**
- `num_filters`：卷积滤波器的数量
- `kernel_size`：卷积核大小
- `num_layers`：转换层数
- `dropout`：辍学率

**用途：**
<<<代码块_2>>>

**时间卷积网络** (`TCN`)
- 用于长程依赖的扩张卷积
- 因果卷积（没有未来信息泄漏）
- 对于长序列有效
- 适合时间序列预测

**优点：**
- 捕获长期依赖关系
- 可并行化（比 RNN 更快）
- 稳定的梯度

### 循环神经网络

**RNN** (`RNN`)
- 基本的循环架构
- 支持 LSTM、GRU、RNN 变体
- 顺序处理
- 捕获时间依赖性

**参数：**
- `rnn_type`：“LSTM”、“GRU”或“RNN”
- `hidden_dim`：隐藏状态维度
- `num_layers`：循环层数
- `dropout`：辍学率
- `bidirectional`：使用双向 RNN

**用途：**
<<<代码块_3>>>

**最适合：**
- 连续的临床事件
- 时间模式学习
- 可变长度序列

### 变压器模型

**变压器** (`Transformer`)
- 自注意力机制
- 序列的并行处理
- 最先进的性能
- 对远程依赖有效

**架构：**
- 多头自注意力
- 位置嵌入
- 前馈网络
- 层标准化

**参数：**
- `num_heads`：注意力头的数量
- `num_layers`：变压器层数
- `hidden_dim`：隐藏维度
- `dropout`：辍学率
- `max_seq_length`：最大序列长度

**用途：**
<<<代码块_4>>>

**变压器模型** (`TransformersModel`)
- 与 HuggingFace 变压器集成
- 用于临床文本的预训练语言模型
- 针对医疗保健任务进行微调
- 示例：BERT、RoBERTa、BioClinicalBERT

**用途：**
<<<代码块_5>>>

### 图神经网络

**GNN** (`GNN`)
- 基于图的学习
- 建模实体之间的关系
- 支持GAT（图注意力）和GCN（图卷积）

**使用案例：**
- 药物间相互作用
- 患者相似网络
- 知识图谱集成
- 合并症关系

**参数：**
- `gnn_type`：“GAT”或“GCN”
- `hidden_dim`：隐藏维度
- `num_layers`：GNN 层数
- `dropout`：辍学率
- `num_heads`：注意头（适用于 GAT）

**用途：**
<<<代码块_6>>>

## 医疗保健特定模型

### 可解释的临床模型

**保留** (`RETAIN`)
- 逆时注意力机制
- 高度可解释的预测
- 访问级别和活动级别的关注
- 识别有影响的临床事件

**主要特点：**
- 两级关注（访问量和功能）
- 时间衰减建模
- 具有临床意义的解释
- 发表于 NeurIPS 2016

**用途：**
```python
from pyhealth.models import RETAIN

model = RETAIN(
    dataset=sample_dataset,
    feature_keys=["diagnoses", "medications"],
    mode="binary",
    hidden_dim=128
)

# Get attention weights for interpretation
outputs = model(batch)
visit_attention = outputs["visit_attention"]
feature_attention = outputs["feature_attention"]
```

**最适合：**
- 死亡率预测
- 再入院预测
- 临床风险评分
- 可解释的预测

**AdaCare** (`AdaCare`)
- 具有特征校准的自适应护理模型
- 针对特定疾病的关注
- 处理不规则的时间间隔
- 可解释的特征重要性

**ConCare** (`ConCare`)
- 交叉访问卷积注意力
- 时间卷积特征提取
- 多层次注意力机制
- 适合纵向 EHR 建模

### 药物推荐模型

**GAMENet** (`GAMENet`)
- 基于图表的药物推荐
- 药物相互作用建模
- 患者病史记忆网络
- 多跳推理

**架构：**
- 药物知识图谱
- 记忆增强神经网络
- DDI感知预测

**用途：**
```python
from pyhealth.models import GAMENet

model = GAMENet(
    dataset=sample_dataset,
    feature_keys=["diagnoses", "medications"],
    mode="multilabel",
    embedding_dim=128,
    ddi_adj_path="/path/to/ddi_adjacency_matrix.pkl"
)
```

**微米** (`MICRON`)
- 具有 DDI 限制的药物推荐
- 交互感知预测
- 注重安全的药物选择

**安全药物** (`SafeDrug`)
- 安全意识药物推荐
- 分子结构整合
- DDI约束优化
- 平衡功效和安全性

**主要特点：**
- 分子图编码
- DDI图神经网络
- 安全强化学习
- 发表于 KDD 2021

**用途：**
```python
from pyhealth.models import SafeDrug

model = SafeDrug(
    dataset=sample_dataset,
    feature_keys=["diagnoses", "medications"],
    mode="multilabel",
    ddi_adj_path="/path/to/ddi_matrix.pkl",
    molecule_path="/path/to/molecule_graphs.pkl"
)
```

**MoleRec** (`MoleRec`)
- 分子水平药物推荐
- 子结构推理
- 细粒度的药物选择

### 疾病进展模型

**StageNet** (`StageNet`)
- 疾病阶段感知预测
- 自动学习临床阶段
- 阶段自适应特征提取
- 有效监测慢性病

**架构：**
- 阶段感知 LSTM
- 动态阶段转换
- 时间衰减机制

**用途：**
```python
from pyhealth.models import StageNet

model = StageNet(
    dataset=sample_dataset,
    feature_keys=["diagnoses", "medications"],
    mode="binary",
    hidden_dim=128,
    num_stages=3,
    chunk_size=128
)
```

**最适合：**
- ICU死亡率预测
- 慢性疾病进展
- 时变风险评估

**更深** (`Deepr`)
- 深度循环架构
- 医学概念嵌入
- 时间模式学习
- 发表于 JAMIA

### 高级序列模型

**代理** (`Agent`)
- 基于强化学习
- 治疗建议
- 行动价值优化
- 顺序决策的策略学习

**掌握** (`GRASP`)
- 基于图的序列模式
- 结构事件关系
- 分层表示学习

**SparcNet** (`SparcNet`)
- 稀疏的临床网络
- 高效的特征选择
- 降低计算成本
- 可解释的预测

**ContraWR** (`ContraWR`)
- 对比学习法
- 自监督预训练
- 稳健的表述
- 有限的标记数据场景

### 医疗实体链接

**MedLink** (`MedLink`)
- 链接到知识库的医疗实体
- 临床概念标准化
- UMLS集成
- 实体消歧

### 生成模型

**GAN** (`GAN`)
- 生成对抗网络
- 综合 EHR 数据生成
- 保护隐私的数据共享
- 针对罕见情况的增强

**VAE** (`VAE`)
- 变分自动编码器
- 患者代表学习
- 异常检测
- 潜在太空探索

### 健康的社会决定因素

**SDOH** (`SDOH`)
- 社会决定因素整合
- 多模态预测
- 解决健康差异
- 结合临床和社会数据

## 型号选择指南

### 按任务类型

**二元分类**（死亡率、再入院）
- 开始于：Logistic 回归（基线）
- 标准：RNN、Transformer
- 可解释：RETAIN、AdaCare
- 高级：StageNet

**多标签分类**（药物推荐）
- 标准：CNN、RNN
- 医疗保健专用：GAMENet、SafeDrug、MICRON、MoleRec
- 基于图：GNN

**回归**（停留时间）
- 开始于：MLP（基线）
- 顺序：RNN、TCN
- 高级：变压器

**多类分类**（医学编码、专业）
- 标准：CNN、RNN、Transformer
- 基于文本：TransformersModel（BERT 变体）

### 按数据类型
**连续事件**（诊断、药物治疗、手术）
- RNN、LSTM、GRU
- 变压器
- 保留、AdaCare、ConCare

**时间序列信号**（脑电图、心电图）
- 美国有线电视新闻网、TCN
-循环神经网络
- 变压器

**文本**（临床笔记）
- TransformersModel（ClinicalBERT、BioBERT）
- CNN 用于较短的文本
- 用于连续文本的 RNN

**图表**（药物相互作用、患者网络）
- GNN（GAT、GCN）
- GAMENet、SafeDrug

**图像**（X 射线、CT 扫描）
- CNN（ResNet、通过 TransformersModel 的 DenseNet）
- 视觉变形金刚

### 通过可解释性需求

**需要高可解释性：**
- 逻辑回归
- 保留
- 阿达护理
- 斯帕克网

**中等可解释性：**
- CNN（过滤器可视化）
- Transformer（注意力可视化）
- GNN（图注意力）

**黑盒可接受：**
- 深度RNN模型
- 复杂的合奏

## 培训注意事项

### 超参数调整

**嵌入尺寸：**
- 小数据集：64-128
- 大型数据集：128-256
- 复杂任务：256-512

**隐藏维度：**
- 与 embedding_dim 成比例
- 通常为 1-2x embedding_dim

**层数：**
- 从2-3层开始
- 更深入的复杂图案
- 注意过度拟合

**退学：**
- 从 0.5 开始
- 如果欠拟合则减少 (0.1-0.3)
- 如果过度拟合则增加 (0.5-0.7)

### 计算要求

**内存（GPU）：**
- CNN：低到中等
- RNN：中等（取决于序列长度）
- 变压器：高（序列长度的二次方）
- GNN：中到高（取决于图大小）

**训练速度：**
- 最快：Logistic 回归、MLP、CNN
- 中等：RNN、GNN
- 较慢：Transformer（但可并行）

### 最佳实践

1. **从简单的基线开始**（逻辑回归，MLP）
2. **根据数据可用性使用适当的功能键**
3. **将模式与任务输出匹配**（二元、多类、多标签、回归）
4. **考虑临床部署的可解释性要求**
5. **验证保留的测试集**以获得真实的性能
6. **监控过度拟合**，尤其是复杂模型
7. **尽可能使用预训练模型** (TransformersModel)
8. **考虑部署的计算限制**

## 工作流程示例

```python
from pyhealth.datasets import MIMIC4Dataset
from pyhealth.tasks import mortality_prediction_mimic4_fn
from pyhealth.models import Transformer
from pyhealth.trainer import Trainer

# 1. Prepare data
dataset = MIMIC4Dataset(root="/path/to/data")
sample_dataset = dataset.set_task(mortality_prediction_mimic4_fn)

# 2. Initialize model
model = Transformer(
    dataset=sample_dataset,
    feature_keys=["diagnoses", "medications", "procedures"],
    mode="binary",
    embedding_dim=128,
    num_heads=8,
    num_layers=3,
    dropout=0.3
)

# 3. Train model
trainer = Trainer(model=model)
trainer.train(
    train_dataloader=train_loader,
    val_dataloader=val_loader,
    epochs=50,
    monitor="pr_auc_score",
    monitor_criterion="max"
)

# 4. Evaluate
results = trainer.evaluate(test_loader)
print(results)
```