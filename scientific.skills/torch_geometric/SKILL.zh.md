<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 火炬几何
描述：“图神经网络 (PyG)。节点/图分类、链接预测、GCN、GAT、GraphSAGE、异构图、分子属性预测，用于几何深度学习。”
---

# PyTorch 几何 (PyG)

## 概述

PyTorch Geometric 是一个基于 PyTorch 构建的库，用于开发和训练图神经网络 (GNN)。将此技能应用于图形和不规则结构的深度学习，包括小批量处理、多 GPU 训练和几何深度学习应用。

## 何时使用此技能

在处理以下情况时应使用此技能：
- **基于图的机器学习**：节点分类、图分类、链接预测
- **分子性质预测**：药物发现、化学性质预测
- **社交网络分析**：社区检测、影响力预测
- **引文网络**：论文分类、推荐系统
- **3D 几何数据**：点云、网格、分子结构
- **异构图**：多类型节点和边（例如知识图）
- **大规模图学习**：邻居采样、分布式训练

## 快速入门

### 安装

```bash
uv pip install torch_geometric
```

对于其他依赖项（稀疏操作、集群）：
<<<代码块_1>>>

### 基本图形创建

<<<代码块_2>>>

### 加载基准数据集

<<<代码块_3>>>

## 核心概念

### 数据结构

PyG 使用具有以下关键属性的 `torch_geometric.data.Data` 类表示图形：

- **`data.x`**：节点特征矩阵`[num_nodes, num_node_features]`
- **`data.edge_index`**：COO 格式的图形连接`[2, num_edges]`
- **`data.edge_attr`**：边缘特征矩阵`[num_edges, num_edge_features]`（可选）
- **`data.y`**：节点或图形的目标标签
- **`data.pos`**：节点空间位置`[num_nodes, num_dimensions]`（可选）
- **自定义属性**：可以添加任何属性（例如，`data.train_mask`、`data.batch`）

**重要**：这些属性不是强制性的 - 根据需要使用自定义属性扩展数据对象。

### 边缘索引格式

边以 COO（坐标）格式存储为 `[2, num_edges]` 张量：
- 第一行：源节点索引
- 第二行：目标节点索引

<<<代码块_4>>>

### 小批量处理

PyG 通过创建块对角邻接矩阵、将多个图连接成一个大的断开连接图来处理批处理：

- 邻接矩阵对角堆叠
- 节点特征沿节点维度串联
- `batch` 向量将每个节点映射到其源图
- 无需填充——计算效率高

<<<代码块_5>>>

## 构建图神经网络

### 消息传递范式

PyG 中的 GNN 遵循邻域聚合方案：
1. 变换节点特征
2. 沿边缘传播消息
3. 聚合来自邻居的消息
4. 更新节点表示

### 使用预构建层

PyG 提供 40 多个卷积层。常见的包括：

**GCNConv**（图卷积网络）：
<<<代码块_6>>>

**GATConv**（图注意力网络）：
```python
from torch_geometric.nn import GATConv

class GAT(torch.nn.Module):
    def __init__(self, num_features, num_classes):
        super().__init__()
        self.conv1 = GATConv(num_features, 8, heads=8, dropout=0.6)
        self.conv2 = GATConv(8 * 8, num_classes, heads=1, concat=False, dropout=0.6)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = F.dropout(x, p=0.6, training=self.training)
        x = F.elu(self.conv1(x, edge_index))
        x = F.dropout(x, p=0.6, training=self.training)
        x = self.conv2(x, edge_index)
        return F.log_softmax(x, dim=1)
```

**图圣人**：
```python
from torch_geometric.nn import SAGEConv

class GraphSAGE(torch.nn.Module):
    def __init__(self, num_features, num_classes):
        super().__init__()
        self.conv1 = SAGEConv(num_features, 64)
        self.conv2 = SAGEConv(64, num_classes)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        return F.log_softmax(x, dim=1)
```

### 自定义消息传递层

对于自定义图层，继承自 `MessagePassing`：

```python
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops, degree

class CustomConv(MessagePassing):
    def __init__(self, in_channels, out_channels):
        super().__init__(aggr='add')  # "add", "mean", or "max"
        self.lin = torch.nn.Linear(in_channels, out_channels)

    def forward(self, x, edge_index):
        # Add self-loops to adjacency matrix
        edge_index, _ = add_self_loops(edge_index, num_nodes=x.size(0))

        # Transform node features
        x = self.lin(x)

        # Compute normalization
        row, col = edge_index
        deg = degree(col, x.size(0), dtype=x.dtype)
        deg_inv_sqrt = deg.pow(-0.5)
        norm = deg_inv_sqrt[row] * deg_inv_sqrt[col]

        # Propagate messages
        return self.propagate(edge_index, x=x, norm=norm)

    def message(self, x_j, norm):
        # x_j: features of source nodes
        return norm.view(-1, 1) * x_j
```

关键方法：
- **`forward()`**：主入口点
- **`message()`**：构造从源节点到目标节点的消息
- **`aggregate()`**：聚合消息（通常不覆盖 - 设置`aggr`参数）
- **`update()`**：聚合后更新节点嵌入

**变量命名约定**：将 `_i` 或 `_j` 附加到张量名称会自动将它们映射到目标或源节点。

## 使用数据集

### 加载内置数据集

PyG 提供了广泛的基准数据集：

```python
# Citation networks (node classification)
from torch_geometric.datasets import Planetoid
dataset = Planetoid(root='/tmp/Cora', name='Cora')  # or 'CiteSeer', 'PubMed'

# Graph classification
from torch_geometric.datasets import TUDataset
dataset = TUDataset(root='/tmp/ENZYMES', name='ENZYMES')

# Molecular datasets
from torch_geometric.datasets import QM9
dataset = QM9(root='/tmp/QM9')

# Large-scale datasets
from torch_geometric.datasets import Reddit
dataset = Reddit(root='/tmp/Reddit')
```

检查 `references/datasets_reference.md` 以获得完整列表。

### 创建自定义数据集

对于适合内存的数据集，继承自 `InMemoryDataset`：

```python
from torch_geometric.data import InMemoryDataset, Data
import torch

class MyOwnDataset(InMemoryDataset):
    def __init__(self, root, transform=None, pre_transform=None):
        super().__init__(root, transform, pre_transform)
        self.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return ['my_data.csv']  # Files needed in raw_dir

    @property
    def processed_file_names(self):
        return ['data.pt']  # Files in processed_dir

    def download(self):
        # Download raw data to self.raw_dir
        pass

    def process(self):
        # Read data, create Data objects
        data_list = []

        # Example: Create a simple graph
        edge_index = torch.tensor([[0, 1], [1, 0]], dtype=torch.long)
        x = torch.randn(2, 16)
        y = torch.tensor([0], dtype=torch.long)

        data = Data(x=x, edge_index=edge_index, y=y)
        data_list.append(data)

        # Apply pre_filter and pre_transform
        if self.pre_filter is not None:
            data_list = [d for d in data_list if self.pre_filter(d)]

        if self.pre_transform is not None:
            data_list = [self.pre_transform(d) for d in data_list]

        # Save processed data
        self.save(data_list, self.processed_paths[0])
```

对于无法容纳在内存中的大型数据集，请继承 `Dataset` 并实现 `len()` 和 `get(idx)`。

### 从 CSV 加载图表

```python
import pandas as pd
import torch
from torch_geometric.data import HeteroData

# Load nodes
nodes_df = pd.read_csv('nodes.csv')
x = torch.tensor(nodes_df[['feat1', 'feat2']].values, dtype=torch.float)

# Load edges
edges_df = pd.read_csv('edges.csv')
edge_index = torch.tensor([edges_df['source'].values,
                           edges_df['target'].values], dtype=torch.long)

data = Data(x=x, edge_index=edge_index)
```

## 培训工作流程

### 节点分类（单图）

```python
import torch
import torch.nn.functional as F
from torch_geometric.datasets import Planetoid

# Load dataset
dataset = Planetoid(root='/tmp/Cora', name='Cora')
data = dataset[0]

# Create model
model = GCN(dataset.num_features, dataset.num_classes)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)

# Training
model.train()
for epoch in range(200):
    optimizer.zero_grad()
    out = model(data)
    loss = F.nll_loss(out[data.train_mask], data.y[data.train_mask])
    loss.backward()
    optimizer.step()

    if epoch % 10 == 0:
        print(f'Epoch {epoch}, Loss: {loss.item():.4f}')

# Evaluation
model.eval()
pred = model(data).argmax(dim=1)
correct = (pred[data.test_mask] == data.y[data.test_mask]).sum()
acc = int(correct) / int(data.test_mask.sum())
print(f'Test Accuracy: {acc:.4f}')
```

### 图分类（多图）

```python
from torch_geometric.datasets import TUDataset
from torch_geometric.loader import DataLoader
from torch_geometric.nn import global_mean_pool

class GraphClassifier(torch.nn.Module):
    def __init__(self, num_features, num_classes):
        super().__init__()
        self.conv1 = GCNConv(num_features, 64)
        self.conv2 = GCNConv(64, 64)
        self.lin = torch.nn.Linear(64, num_classes)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        x = F.relu(x)

        # Global pooling (aggregate node features to graph-level)
        x = global_mean_pool(x, batch)

        x = self.lin(x)
        return F.log_softmax(x, dim=1)

# Load dataset
dataset = TUDataset(root='/tmp/ENZYMES', name='ENZYMES')
loader = DataLoader(dataset, batch_size=32, shuffle=True)

model = GraphClassifier(dataset.num_features, dataset.num_classes)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

# Training
model.train()
for epoch in range(100):
    total_loss = 0
    for batch in loader:
        optimizer.zero_grad()
        out = model(batch)
        loss = F.nll_loss(out, batch.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()

    if epoch % 10 == 0:
        print(f'Epoch {epoch}, Loss: {total_loss / len(loader):.4f}')
```

### 具有邻域采样的大型图

对于大图，使用 `NeighborLoader` 对子图进行采样：
```python
from torch_geometric.loader import NeighborLoader

# Create a neighbor sampler
train_loader = NeighborLoader(
    data,
    num_neighbors=[25, 10],  # Sample 25 neighbors for 1st hop, 10 for 2nd hop
    batch_size=128,
    input_nodes=data.train_mask,
)

# Training
model.train()
for batch in train_loader:
    optimizer.zero_grad()
    out = model(batch)
    # Only compute loss on seed nodes (first batch_size nodes)
    loss = F.nll_loss(out[:batch.batch_size], batch.y[:batch.batch_size])
    loss.backward()
    optimizer.step()
```

**重要**：
- 输出子图是有向的
- 节点索引被重新标记（0 到batch.num_nodes - 1）
- 仅使用种子节点预测进行损失计算
- 超过 2-3 跳的采样通常是不可行的

## 高级功能

### 异构图

对于具有多个节点和边类型的图，请使用 `HeteroData`：

```python
from torch_geometric.data import HeteroData

data = HeteroData()

# Add node features for different types
data['paper'].x = torch.randn(100, 128)  # 100 papers with 128 features
data['author'].x = torch.randn(200, 64)  # 200 authors with 64 features

# Add edges for different types (source_type, edge_type, target_type)
data['author', 'writes', 'paper'].edge_index = torch.randint(0, 200, (2, 500))
data['paper', 'cites', 'paper'].edge_index = torch.randint(0, 100, (2, 300))

print(data)
```

将同类模型转换为异构模型：

```python
from torch_geometric.nn import to_hetero

# Define homogeneous model
model = GNN(...)

# Convert to heterogeneous
model = to_hetero(model, data.metadata(), aggr='sum')

# Use as normal
out = model(data.x_dict, data.edge_index_dict)
```

或者使用 `HeteroConv` 进行自定义边类型特定操作：

```python
from torch_geometric.nn import HeteroConv, GCNConv, SAGEConv

class HeteroGNN(torch.nn.Module):
    def __init__(self, metadata):
        super().__init__()
        self.conv1 = HeteroConv({
            ('paper', 'cites', 'paper'): GCNConv(-1, 64),
            ('author', 'writes', 'paper'): SAGEConv((-1, -1), 64),
        }, aggr='sum')

        self.conv2 = HeteroConv({
            ('paper', 'cites', 'paper'): GCNConv(64, 32),
            ('author', 'writes', 'paper'): SAGEConv((64, 64), 32),
        }, aggr='sum')

    def forward(self, x_dict, edge_index_dict):
        x_dict = self.conv1(x_dict, edge_index_dict)
        x_dict = {key: F.relu(x) for key, x in x_dict.items()}
        x_dict = self.conv2(x_dict, edge_index_dict)
        return x_dict
```

### 变换

应用变换来修改图结构或特征：

```python
from torch_geometric.transforms import NormalizeFeatures, AddSelfLoops, Compose

# Single transform
transform = NormalizeFeatures()
dataset = Planetoid(root='/tmp/Cora', name='Cora', transform=transform)

# Compose multiple transforms
transform = Compose([
    AddSelfLoops(),
    NormalizeFeatures(),
])
dataset = Planetoid(root='/tmp/Cora', name='Cora', transform=transform)
```

常见变换：
- **结构**：`ToUndirected`、`AddSelfLoops`、`RemoveSelfLoops`、`KNNGraph`、`RadiusGraph`
- **功能**：`NormalizeFeatures`、`NormalizeScale`、`Center`
- **采样**：`RandomNodeSplit`、`RandomLinkSplit`
- **位置编码**：`AddLaplacianEigenvectorPE`、`AddRandomWalkPE`

完整列表请参见`references/transforms_reference.md`。

### 模型可解释性

PyG 提供了可解释性工具来理解模型预测：

```python
from torch_geometric.explain import Explainer, GNNExplainer

# Create explainer
explainer = Explainer(
    model=model,
    algorithm=GNNExplainer(epochs=200),
    explanation_type='model',  # or 'phenomenon'
    node_mask_type='attributes',
    edge_mask_type='object',
    model_config=dict(
        mode='multiclass_classification',
        task_level='node',
        return_type='log_probs',
    ),
)

# Generate explanation for a specific node
node_idx = 10
explanation = explainer(data.x, data.edge_index, index=node_idx)

# Visualize
print(f'Node {node_idx} explanation:')
print(f'Important edges: {explanation.edge_mask.topk(5).indices}')
print(f'Important features: {explanation.node_mask[node_idx].topk(5).indices}')
```

### 池化操作

对于层次图表示：

```python
from torch_geometric.nn import TopKPooling, global_mean_pool

class HierarchicalGNN(torch.nn.Module):
    def __init__(self, num_features, num_classes):
        super().__init__()
        self.conv1 = GCNConv(num_features, 64)
        self.pool1 = TopKPooling(64, ratio=0.8)
        self.conv2 = GCNConv(64, 64)
        self.pool2 = TopKPooling(64, ratio=0.8)
        self.lin = torch.nn.Linear(64, num_classes)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        x = F.relu(self.conv1(x, edge_index))
        x, edge_index, _, batch, _, _ = self.pool1(x, edge_index, None, batch)

        x = F.relu(self.conv2(x, edge_index))
        x, edge_index, _, batch, _, _ = self.pool2(x, edge_index, None, batch)

        x = global_mean_pool(x, batch)
        x = self.lin(x)
        return F.log_softmax(x, dim=1)
```

## 常见模式和最佳实践

### 检查图形属性

```python
# Undirected check
from torch_geometric.utils import is_undirected
print(f"Is undirected: {is_undirected(data.edge_index)}")

# Connected components
from torch_geometric.utils import connected_components
print(f"Connected components: {connected_components(data.edge_index)}")

# Contains self-loops
from torch_geometric.utils import contains_self_loops
print(f"Has self-loops: {contains_self_loops(data.edge_index)}")
```

### GPU 训练

```python
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = model.to(device)
data = data.to(device)

# For DataLoader
for batch in loader:
    batch = batch.to(device)
    # Train...
```

### 保存和加载模型

```python
# Save
torch.save(model.state_dict(), 'model.pth')

# Load
model = GCN(num_features, num_classes)
model.load_state_dict(torch.load('model.pth'))
model.eval()
```

### 层功能

选择层时，请考虑以下功能：
- **SparseTensor**：支持高效的稀疏矩阵运算
- **edge_weight**：处理一维边缘权重
- **edge_attr**：处理多维边缘特征
- **二分图**：适用于二分图（不同的源/目标维度）
- **Lazy**：启用初始化而不指定输入维度

请参阅 `references/layer_capabilities.md` 处的 GNN 备忘单。

## 资源

### 捆绑参考资料

该技能包括详细的参考文档：

- **`references/layers_reference.md`**：所有 40 多个 GNN 层的完整列表以及描述和功能
- **`references/datasets_reference.md`**：按类别组织的综合数据集目录
- **`references/transforms_reference.md`**：所有可用的转换及其用例
- **`references/api_patterns.md`**：常见 API 模式和编码示例

### 脚本

`scripts/` 中提供了实用程序脚本：

- **`scripts/visualize_graph.py`**：使用networkx和matplotlib可视化图形结构
- **`scripts/create_gnn_template.py`**：为常见的 GNN 架构生成样板代码
- **`scripts/benchmark_model.py`**：标准数据集上的基准模型性能

直接执行脚本或读取它们以获取实现模式。

### 官方资源

- **文档**：https://pytorch-geometric.readthedocs.io/
- **GitHub**：https://github.com/pyg-team/pytorch_geometric
- **教程**：https://pytorch-geometric.readthedocs.io/en/latest/get_started/introduction.html
- **示例**：https://github.com/pyg-team/pytorch_geometric/tree/master/examples