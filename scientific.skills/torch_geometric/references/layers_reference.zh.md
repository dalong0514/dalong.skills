<!-- 此文件由机器翻译自 layers_reference.md -->

# PyTorch 几何神经网络层参考

本文档提供了 `torch_geometric.nn` 中可用的所有神经网络层的综合参考。

## 层能力标志

选择层时，请考虑以下功能标志：

- **SparseTensor**：支持 `torch_sparse.SparseTensor` 格式以实现高效的稀疏操作
- **edge_weight**：处理一维边缘权重数据
- **edge_attr**：处理多维边缘特征信息
- **二分**：适用于二分图（不同的源/目标节点维度）
- **静态**：在具有批处理节点功能的静态图上运行
- **Lazy**：启用初始化而不指定输入通道尺寸

## 卷积层

### 标准图卷积

**GCNConv** - 图卷积网络层
- 通过对称归一化实现谱图卷积
- 支持：SparseTensor、edge_weight、Bipartite、Lazy
- 用于：引文网络、社交网络、一般图学习
- 示例：`GCNConv(in_channels, out_channels, improved=False, cached=True)`

**SAGEConv** - GraphSAGE 层
- 通过邻域采样和聚合进行归纳学习
- 支持：SparseTensor、Bipartite、Lazy
- 用于：大图、归纳学习、异构特征
- 示例：`SAGEConv(in_channels, out_channels, aggr='mean')`

**GATConv** - 图注意力网络层
- 用于自适应邻居加权的多头注意机制
- 支持：SparseTensor、edge_attr、Bipartite、Static、Lazy
- 用于：需要可变邻居重要性的任务
- 示例：`GATConv(in_channels, out_channels, heads=8, dropout=0.6)`

**GraphConv** - 简单图卷积（Morris 等人）
- 具有可选边缘权重的基本消息传递
- 支持：SparseTensor、edge_weight、Bipartite、Lazy
- 用于：基线模型、简单的图形结构
- 示例：`GraphConv(in_channels, out_channels, aggr='add')`

**GINConv** - 图同构网络层
- 用于图同构测试的最强大的 GNN
- 支持：双方
- 用于：图分类、分子特性预测
- 示例：`GINConv(nn.Sequential(nn.Linear(in_channels, out_channels), nn.ReLU()))`

**TransformerConv** - Graph Transformer 层
- 将图结构与变压器注意力相结合
- 支持：SparseTensor、Bipartite、Lazy
- 用于：远程依赖、复杂图表
- 示例：`TransformerConv(in_channels, out_channels, heads=8, beta=True)`

**ChebConv** - 切比雪夫谱图卷积
- 使用切比雪夫多项式进行有效的光谱过滤
- 支持：SparseTensor、edge_weight、Bipartite、Lazy
- 用于：谱图学习、高效卷积
- 示例：`ChebConv(in_channels, out_channels, K=3)`

**SGConv** - 简化图卷积
- 预先计算固定数量的传播步骤
- 支持：SparseTensor、edge_weight、Bipartite、Lazy
- 用于：快速训练、浅层模型
- 示例：`SGConv(in_channels, out_channels, K=2)`

**APPNP** - 神经预测的近似个性化传播
- 将特征转换与传播分开
- 支持：SparseTensor、edge_weight、Lazy
- 用于：深度传播而不过度平滑
- 示例：`APPNP(K=10, alpha=0.1)`

**ARMAConv** - ARMA 图卷积
- 使用 ARMA 过滤器进行图形过滤
- 支持：SparseTensor、edge_weight、Bipartite、Lazy
- 用于：高级光谱方法
- 示例：`ARMAConv(in_channels, out_channels, num_stacks=3, num_layers=2)`

**GATv2Conv** - 改进的图注意力网络
- 修复了 GAT 中的静态注意力计算问题
- 支持：SparseTensor、edge_attr、Bipartite、Static、Lazy
- 用于：比原始 GAT 更好的注意力学习
- 示例：`GATv2Conv(in_channels, out_channels, heads=8)`

**SuperGATConv** - 自监督图注意力
- 增加自我监督注意力机制
- 支持：SparseTensor、edge_attr、Bipartite、Static、Lazy
- 用于：自我监督学习、有限标签
- 示例：`SuperGATConv(in_channels, out_channels, heads=8)`

**GMMConv** - 高斯混合模型卷积
- 在伪坐标空间中使用高斯核
- 支持：双方
- 用于：点云、空间数据
- 示例：`GMMConv(in_channels, out_channels, dim=3, kernel_size=5)`

**SplineConv** - 基于样条的卷积
- 用于空间过滤的 B 样条基函数
- 支持：双方
- 用于：不规则网格、连续空间
- 示例：`SplineConv(in_channels, out_channels, dim=2, kernel_size=5)`

**NNConv** - 神经网络卷积
- 神经网络处理的边缘特征
- 支持：edge_attr、二分
- 用于：丰富的边缘特征、分子图
- 示例：`NNConv(in_channels, out_channels, nn=edge_nn, aggr='mean')`

**CGConv** - 水晶图卷积
- 专为结晶材料而设计
- 支持：双方
- 用于：材料科学、晶体结构
- 示例：`CGConv(in_channels, dim=3, batch_norm=True)`

**EdgeConv** - 边缘卷积（动态图 CNN）
- 基于特征空间动态计算边缘
支持：静态
- 用于：点云、动态图
- 示例：`EdgeConv(nn=edge_nn, aggr='max')`

**PointNetConv** - PointNet++ 卷积
- 点云的局部和全局特征学习
- 用于：3D 点云处理
- 示例：`PointNetConv(local_nn, global_nn)`

**ResGatedGraphConv** - 残差门控图卷积
- 具有残余连接的门控机制
- 支持：edge_attr、Bipartite、Lazy
- 用于：深度 GNN、复杂特征
- 示例：`ResGatedGraphConv(in_channels, out_channels)`

**GENConv** - 广义图卷积
- 概括多个 GNN 变体
- 支持：SparseTensor、edge_weight、edge_attr、Bipartite、Lazy
- 用于：灵活的架构探索
- 示例：`GENConv(in_channels, out_channels, aggr='softmax', num_layers=2)`

**FiLMConv** - 特征线性调制
- 全局特征的条件
- 支持：双向、惰性
- 用于：条件生成、多任务学习
- 示例：`FiLMConv(in_channels, out_channels, num_relations=5)`

**PANConv** - 路径注意力网络
- 注意多跳路径
- 支持：SparseTensor、Lazy
- 用于：复杂的连接模式
- 示例：`PANConv(in_channels, out_channels, filter_size=3)`

**ClusterGCNConv** - Cluster-GCN 卷积
- 通过图聚类进行高效训练
- 支持：edge_attr、Lazy
- 用于：非常大的图表
- 示例：`ClusterGCNConv(in_channels, out_channels)`

**MFConv** - 多尺度特征卷积
- 聚合多个尺度的特征
- 支持：SparseTensor、Lazy
- 用于：多尺度图案
- 示例：`MFConv(in_channels, out_channels)`

**RGCNConv** - 关系图卷积
- 处理多种边缘类型
- 支持：SparseTensor、edge_weight、Lazy
- 用于：知识图、异构图
- 示例：`RGCNConv(in_channels, out_channels, num_relations=10)`

**FAConv** - 频率自适应卷积
- 谱域自适应滤波
- 支持：SparseTensor、Lazy
- 用于：光谱图学习
- 示例：`FAConv(in_channels, eps=0.1, dropout=0.5)`

### 分子和 3D 卷积

**SchNet** - 连续滤波器卷积层
- 专为分子动力学而设计
- 用于：分子特性预测、3D 分子
- 示例：`SchNet(hidden_channels=128, num_filters=64, num_interactions=6)`

**DimeNet** - 定向消息传递
- 使用方向信息和角度
- 用于：3D 分子结构、化学性质
- 示例：`DimeNet(hidden_channels=128, out_channels=1, num_blocks=6)`

**PointTransformerConv** - 点云转换器
- 3D 点云转换器
- 用于：3D 视觉、点云分割
- 示例：`PointTransformerConv(in_channels, out_channels)`

### 超图卷积

**HypergraphConv** - 超图卷积
- 在超边（连接多个节点的边）上运行
- 支持：懒惰
- 用于：多向关系、化学反应
- 示例：`HypergraphConv(in_channels, out_channels)`

**HGTConv** - 异构图转换器
- 多种类型异构图的转换器
- 支持：懒惰
- 用于：异构网络、知识图
- 示例：`HGTConv(in_channels, out_channels, metadata, heads=8)`

## 聚合运算符

**Aggr** - 基础聚合类
- 跨节点灵活聚合

**SumAggregation** - 求和聚合
- 示例：`SumAggregation()`

**MeanAggregation** - 平均聚合
- 示例：`MeanAggregation()`

**MaxAggregation** - 最大聚合
- 示例：`MaxAggregation()`

**SoftmaxAggregation** - Softmax 加权聚合
- 可学习的注意力权重
- 示例：`SoftmaxAggregation(learn=True)`

**PowerMeanAggregation** - 幂均值聚合
- 可学习的功率参数
- 示例：`PowerMeanAggregation(learn=True)`

**LSTMAggregation** - 基于 LSTM 的聚合
- 邻居的顺序处理
- 示例：`LSTMAggregation(in_channels, out_channels)`

**SetTransformerAggregation** - 设置 Transformer 聚合
- 用于排列不变聚合的变压器
- 示例：`SetTransformerAggregation(in_channels, out_channels)`

**多重聚合** - 多重聚合
- 结合多种聚合方法
- 示例：`MultiAggregation(['mean', 'max', 'std'])`

## 池化层

### 全球汇集

**global_mean_pool** - 全局均值池
- 平均每个图的节点特征
- 示例：`global_mean_pool(x, batch)`

**global_max_pool** - 全局最大池化
- 每个图的最大节点特征
- 示例：`global_max_pool(x, batch)`

**global_add_pool** - 全局总和池
- 对每个图的节点特征求和
- 示例：`global_add_pool(x, batch)`

**global_sort_pool** - 全局排序池
- 对前 k 个节点进行排序和连接
- 示例：`global_sort_pool(x, batch, k=30)`

**GlobalAttention** - 全球注意力池
- 用于聚合的可学习注意力权重
- 示例：`GlobalAttention(gate_nn)`

**Set2Set** - Set2Set 池化
- 基于LSTM的注意力机制
- 示例：`Set2Set(in_channels, processing_steps=3)`

### 分层池化
**TopKPooling** - Top-k 池化
- 根据投影分数保留前 k 个节点
- 示例：`TopKPooling(in_channels, ratio=0.5)`

**SAGPooling** - 自注意力图池
- 使用自注意力进行节点选择
- 示例：`SAGPooling(in_channels, ratio=0.5)`

**ASAPooling** - 自适应结构感知池
- 结构感知节点选择
- 示例：`ASAPooling(in_channels, ratio=0.5)`

**PANPooling** - 路径注意力池
- 注意池化路径
- 示例：`PANPooling(in_channels, ratio=0.5)`

**EdgePooling** - 边缘收缩池
- 收缩边池
- 示例：`EdgePooling(in_channels)`

**MemPooling** - 基于内存的池
- 可学习的集群分配
- 示例：`MemPooling(in_channels, out_channels, heads=4, num_clusters=10)`

**avg_pool** / **max_pool** - 带聚类的平均/最大池
- 集群内的池节点
- 示例：`avg_pool(cluster, data)`

## 标准化层

**BatchNorm** - 批量归一化
- 跨批次标准化特征
- 示例：`BatchNorm(in_channels)`

**LayerNorm** - 层标准化
- 标准化每个样本的特征
- 示例：`LayerNorm(in_channels)`

**InstanceNorm** - 实例规范化
- 对每个样本和图表进行标准化
- 示例：`InstanceNorm(in_channels)`

**GraphNorm** - 图标准化
- 特定于图的标准化
- 示例：`GraphNorm(in_channels)`

**PairNorm** - 配对标准化
- 防止深度 GNN 中的过度平滑
- 示例：`PairNorm(scale_individually=False)`

**MessageNorm** - 消息规范化
- 在传递过程中规范化消息
- 示例：`MessageNorm(learn_scale=True)`

**DiffGroupNorm** - 可微分组标准化
- 可学习的标准化分组
- 示例：`DiffGroupNorm(in_channels, groups=10)`

## 模型架构

### 预建模型

**GCN** - 完整的图卷积网络
- 带dropout的多层GCN
- 示例：`GCN(in_channels, hidden_channels, num_layers, out_channels)`

**GraphSAGE** - 完整的 GraphSAGE 模型
- 带辍学的多层SAGE
- 示例：`GraphSAGE(in_channels, hidden_channels, num_layers, out_channels)`

**GIN** - 完整的图同构网络
- 用于图分类的多层GIN
- 示例：`GIN(in_channels, hidden_channels, num_layers, out_channels)`

**GAT** - 完整的图注意力网络
- 具有注意力的多层GAT
- 示例：`GAT(in_channels, hidden_channels, num_layers, out_channels, heads=8)`

**PNA** - 主要邻域聚合
- 结合多个聚合器和缩放器
- 示例：`PNA(in_channels, hidden_channels, num_layers, out_channels)`

**EdgeCNN** - 边缘卷积 CNN
- 点云的动态图 CNN
- 示例：`EdgeCNN(out_channels, num_layers=3, k=20)`

### 自动编码器

**GAE** - 图形自动编码器
- 将图编码到潜在空间中
- 示例：`GAE(encoder)`

**VGAE** - 变分图自动编码器
- 概率图编码
- 示例：`VGAE(encoder)`

**ARGA** - 对抗性正则化图自动编码器
- 具有对抗性正则化的 GAE
- 示例：`ARGA(encoder, discriminator)`

**ARGVA** - 对抗正则化变分图自动编码器
- 具有对抗性正则化的 VGAE
- 示例：`ARGVA(encoder, discriminator)`

### 知识图嵌入

**TransE** - 翻译嵌入
- 学习实体和关系嵌入
- 示例：`TransE(num_nodes, num_relations, hidden_channels)`

**RotatE** - 旋转嵌入
- 复杂空间中的嵌入
- 示例：`RotatE(num_nodes, num_relations, hidden_channels)`

**ComplEx** - 复杂嵌入
- 复值嵌入
- 示例：`ComplEx(num_nodes, num_relations, hidden_channels)`

**DistMult** - 双线性对角线模型
- 简化的双线性模型
- 示例：`DistMult(num_nodes, num_relations, hidden_channels)`

## 实用层

**顺序** - 顺序容器
- 链接多层
- 示例：`Sequential('x, edge_index', [(GCNConv(16, 64), 'x, edge_index -> x'), nn.ReLU()])`

**JumpingKnowledge** - 跳跃知识连接
- 组合所有层的表示
- 模式：“cat”、“max”、“lstm”
- 示例：`JumpingKnowledge(mode='cat')`

**DeepGCNLayer** - 深层 GCN 层包装器
- 通过跳过连接启用非常深的 GNN
- 示例：`DeepGCNLayer(conv, norm, act, block='res+', dropout=0.1)`

**MLP** - 多层感知器
- 标准前馈网络
- 示例：`MLP([in_channels, 64, 64, out_channels], dropout=0.5)`

**线性** - 惰性线性层
- 具有延迟初始化的线性变换
- 示例：`Linear(in_channels, out_channels, bias=True)`

## 密集层

对于密集（非稀疏）图表示：

**DenseGCNConv** - 密集 GCN 层
**DenseSAGEConv** - 密集 SAGE 层
**DenseGINConv** - 密集 GIN 层
**DenseGraphConv** - 密集图卷积

这些在处理小型、完全连接或密集表示的图时非常有用。

## 使用提示

1. **从简单开始**：对于大多数任务，从 GCNConv 或 GATConv 开始
2. **考虑数据类型**：将分子层（SchNet、DimeNet）用于 3D 结构
3. **检查功能**：将图层功能与您的数据进行匹配（边缘特征、二分等）
4. **深度网络**：对深度 GNN 使用归一化（PairNorm、LayerNorm）和 JumpingKnowledge
5. **大图**：使用可扩展层（SAGE、Cluster-GCN）和邻居采样
6. **异构**：使用RGCNConv、HGTConv或to_hetero()转换
7. **延迟初始化**：当输入维度变化或未知时使用延迟层

## 常见模式

### 基本 GNN
```python
from torch_geometric.nn import GCNConv, global_mean_pool

class GNN(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, out_channels)

    def forward(self, x, edge_index, batch):
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index)
        return global_mean_pool(x, batch)
```

### 具有归一化的深度 GNN
<<<代码块_1>>>