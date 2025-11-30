<!-- 此文件由机器翻译自 networks.md -->

# 深度学习网络

Aeon 提供专门为时间序列任务设计的神经网络架构。这些网络充当分类、回归、聚类和预测的构建块。

## 核心网络架构

### 卷积网络

**FCNNetwork** - 全卷积网络
- 具有批量归一化功能的三个卷积块
- 用于降维的全局平均池化
- **使用时**：需要简单而有效的 CNN 基线

**ResNetNetwork** - 残差网络
- 具有跳过连接的剩余块
- 防止深层网络中梯度消失
- **使用时**：需要深度网络，训练稳定性很重要

**InceptionNetwork** - Inception 模块
- 使用并行卷积进行多尺度特征提取
- 不同的内核大小捕获不同尺度的模式
- **使用时**：模式存在于多个时间尺度

**TimeCNNNetwork** - 标准 CNN
- 基本卷积架构
- **使用时**：简单的 CNN 就足够了，可解释性很重要

**不相交CNN网络** - 单独的路径
- 不相交的卷积路径
- **使用时**：需要不同的特征提取策略

**DCNNNetwork** - 扩张的 CNN
- 大感受野的扩张卷积
- **使用时**：没有很多层的远程依赖

### 循环网络

**循环网络** - RNN/LSTM/GRU
- 可配置的单元类型（RNN、LSTM、GRU）
- 时间依赖性的顺序建模
- **使用时**：顺序依赖性关键、可变长度系列

### 时间卷积网络

**TCNNetwork** - 时间卷积网络
- 扩张因果卷积
- 感受野大，无复发
- **使用时间**：长序列，需要可并行架构

### 多层感知器

**MLPNetwork** - 基本前馈
- 简单的全连接层
- 处理前压平时间序列
- **使用时间**：需要基线、计算限制或简单模式

## 基于编码器的架构

专为表示学习和聚类而设计的网络。

### 自动编码器变体

**EncoderNetwork** - 通用编码器
- 灵活的编码器结构
- **使用时**：需要自定义编码

**AEFCNNetwork** - 基于 FCN 的自动编码器
- 全卷积编码器-解码器
- **使用时**：需要卷积表示学习

**AEResNetNetwork** - ResNet 自动编码器
- 编码器-解码器中的剩余块
- **使用时间**：具有跳过连接的深度自动编码

**AEDCNNNetwork** - 扩张的 CNN 自动编码器
- 用于压缩的膨胀卷积
- **使用时**：自动编码器中需要大的感受野

**AEDRNNNetwork** - 扩张的 RNN 自动编码器
- 扩张的循环连接
- **使用时间**：具有远程依赖性的顺序模式

**AEBiGRUNetwork** - 双向 GRU
- 双向循环编码
- **使用时间**：来自两个方向的上下文都有帮助

**AEAttentionBiGRUNetwork** - 注意力 + BiGRU
- BiGRU 输出的注意力机制
- **使用时间**：需要关注重要的时间步骤

## 专业架构

**LITENetwork** - 轻量级 Inception Time Ensemble
- 高效的基于初始的架构
- 多元系列的 LITEMV 变体
- **使用时**：需要高效且性能强大

**DeepARNetwork** - 概率预测
- 用于预测的自回归 RNN
- 产生概率预测
- **使用时**：需要预测不确定性量化

## 与估算器一起使用

网络通常在估计器中使用，而不是直接使用：

```python
from aeon.classification.deep_learning import FCNClassifier
from aeon.regression.deep_learning import ResNetRegressor
from aeon.clustering.deep_learning import AEFCNClusterer

# Classification with FCN
clf = FCNClassifier(n_epochs=100, batch_size=16)
clf.fit(X_train, y_train)

# Regression with ResNet
reg = ResNetRegressor(n_epochs=100)
reg.fit(X_train, y_train)

# Clustering with autoencoder
clusterer = AEFCNClusterer(n_clusters=3, n_epochs=100)
labels = clusterer.fit_predict(X_train)
```

## 自定义网络配置

许多网络接受配置参数：

<<<代码块_1>>>

## 基类

- `BaseDeepLearningNetwork` - 所有网络的抽象基础
- `BaseDeepRegressor` - 深度回归的基础
- `BaseDeepClassifier` - 深度分类的基础
- `BaseDeepForecaster` - 深度预测的基础

扩展它们以实现自定义架构。

## 培训注意事项

### 超参数

要调整的关键超参数：

- `n_epochs` - 训练迭代（典型值 50-200）
- `batch_size` - 每批次样本数（典型值 16-64）
- `learning_rate` - 步长 (0.0001-0.01)
- 网络特定：层、过滤器、内核大小

### 回调

许多网络支持训练监控回调：

<<<代码块_2>>>

### GPU 加速

深度学习网络受益于 GPU：

<<<代码块_3>>>

## 架构选择

### 按任务：

**分类**：InceptionNetwork、ResNetNetwork、FCNNetwork
**回归**：InceptionNetwork、ResNetNetwork、TCNNetwork
**预测**：TCNNetwork、DeepARNetwork、RecurrentNetwork
**聚类**：AEFCNNetwork、AEResNetNetwork、AEAttentionBiGRUNetwork

### 按数据特征：

**长序列**：TCNNetwork、DCNNNetwork（扩张卷积）
**短序列**：MLPNetwork、FCNNetwork
**多元**：InceptionNetwork、FCNNetwork、LITENetwork
**可变长度**：带掩码的RecurrentNetwork
**多尺度模式**：InceptionNetwork

### 按计算资源：

**有限计算**：MLPNetwork、LITENetwork
**中等计算**：FCNNetwork、TimeCNNNetwork
**可用的高计算能力**：InceptionNetwork、ResNetNetwork
**GPU可用**：任何深度网络（主要加速）

## 最佳实践

### 1. 数据准备

标准化输入数据：

<<<代码块_4>>>

### 2. 训练/验证分割

使用验证集提前停止：

<<<代码块_5>>>

### 3.从简单开始

先从简单的架构开始，然后再从复杂的架构开始：

1.先尝试MLPNetwork或FCNNetwork
2.如果不够，尝试ResNetNetwork或InceptionNetwork
3. 如果单一模型不足，请考虑集成

### 4. 超参数调优

使用网格搜索或随机搜索：

<<<代码块_6>>>

### 5.正则化

防止过度拟合：
- 使用dropout（如果网络支持）
- 提前停止
- 数据增强（如果可用）
- 降低模型复杂度

### 6. 再现性

设置随机种子：

```python
import numpy as np
import random
import tensorflow as tf

seed = 42
np.random.seed(seed)
random.seed(seed)
tf.random.set_seed(seed)
```