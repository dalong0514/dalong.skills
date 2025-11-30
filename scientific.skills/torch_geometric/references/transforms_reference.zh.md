<!-- 此文件由机器翻译自 transforms_reference.md -->

# PyTorch 几何变换参考

本文档提供了 `torch_geometric.transforms` 中可用的所有转换的综合参考。

## 概述

转换在训练之前或训练期间修改 `Data` 或 `HeteroData` 对象。通过以下方式应用它们：

```python
# During dataset loading
dataset = MyDataset(root='/tmp', transform=MyTransform())

# Apply to individual data
transform = MyTransform()
data = transform(data)

# Compose multiple transforms
from torch_geometric.transforms import Compose
transform = Compose([Transform1(), Transform2(), Transform3()])
```

## 一般变换

### 标准化特征
**目的**：行归一化节点特征，使其总和为 1
**用例**：特征缩放、类概率特征
<<<代码块_1>>>

### 到设备
**用途**：将数据传输到指定设备（CPU/GPU）
**用例**：GPU 训练、设备管理
<<<代码块_2>>>

### 随机节点分割
**目的**：创建训练/验证/测试节点掩码
**用例**：节点分类分割
**参数**：`split='train_rest'`、`num_splits`、`num_val`、`num_test`
<<<代码块_3>>>

### 随机链接分割
**目的**：创建训练/验证/测试边缘分割
**用例**：链接预测
**参数**：`num_val`、`num_test`、`is_undirected`、`split_labels`
<<<代码块_4>>>

### 索引掩码
**用途**：将索引转换为布尔掩码
**用例**：数据预处理
<<<代码块_5>>>

### 掩码到索引
**用途**：将布尔掩码转换为索引
**用例**：数据预处理
<<<代码块_6>>>

### 定点
**目的**：采样固定数量的点
**用例**：点云子采样
**参数**：`num`、`replace`、`allow_duplicates`
```python
from torch_geometric.transforms import FixedPoints
transform = FixedPoints(1024)
```

### 致密
**目的**：转换为稠密邻接矩阵
**用例**：小图，密集操作
```python
from torch_geometric.transforms import ToDense
transform = ToDense(num_nodes=100)
```

### ToSparseTensor
**目的**：将edge_index转换为SparseTensor
**用例**：高效的稀疏操作
**参数**：`remove_edge_index`、`fill_cache`
```python
from torch_geometric.transforms import ToSparseTensor
transform = ToSparseTensor()
```

## 图结构变换

### 至无向
**目的**：将有向图转换为无向图
**用例**：无向图算法
**参数**：`reduce='add'`（如何处理重复的边）
```python
from torch_geometric.transforms import ToUndirected
transform = ToUndirected()
```

### 添加SelfLoops
**目的**：为所有节点添加自循环
**用例**：GCN 式卷积
**参数**：`fill_value`（自循环的边缘属性）
```python
from torch_geometric.transforms import AddSelfLoops
transform = AddSelfLoops()
```

### 删除SelfLoops
**目的**：删除所有自环
**用例**：清理图结构
```python
from torch_geometric.transforms import RemoveSelfLoops
transform = RemoveSelfLoops()
```

### 删除IsolatedNodes
**目的**：删除没有边的节点
**用例**：图形清理
```python
from torch_geometric.transforms import RemoveIsolatedNodes
transform = RemoveIsolatedNodes()
```

### 删除重复边缘
**目的**：删除重复的边
**用例**：图形清理
```python
from torch_geometric.transforms import RemoveDuplicatedEdges
transform = RemoveDuplicatedEdges()
```

### 最大连接组件
**目的**：仅保留最大的连通分量
**用例**：关注主图结构
**参数**：`num_components`（要保留多少个组件）
```python
from torch_geometric.transforms import LargestConnectedComponents
transform = LargestConnectedComponents(num_components=1)
```

### KNNGraph
**目的**：基于k近邻创建边
**用例**：点云、空间数据
**参数**：`k`、`loop`、`force_undirected`、`flow`
```python
from torch_geometric.transforms import KNNGraph
transform = KNNGraph(k=6)
```

### 半径图
**用途**：在半径内创建边缘
**用例**：点云、空间数据
**参数**：`r`、`loop`、`max_num_neighbors`、`flow`
```python
from torch_geometric.transforms import RadiusGraph
transform = RadiusGraph(r=0.1)
```

### 德劳内
**用途**：计算 Delaunay 三角剖分
**用例**：2D/3D 空间图
```python
from torch_geometric.transforms import Delaunay
transform = Delaunay()
```

### 面到边
**用途**：将网格面转换为边
**用例**：网格处理
```python
from torch_geometric.transforms import FaceToEdge
transform = FaceToEdge()
```

### 折线图
**目的**：将图形转换为折线图
**用例**：以边缘为中心的分析
**参数**：`force_directed`
```python
from torch_geometric.transforms import LineGraph
transform = LineGraph()
```

### 游戏开发者大会
**目的**：图扩散卷积预处理
**用例**：改进的消息传递
**参数**：`self_loop_weight`、`normalization_in`、`normalization_out`、`diffusion_kwargs`
```python
from torch_geometric.transforms import GDC
transform = GDC(self_loop_weight=1, normalization_in='sym',
                diffusion_kwargs=dict(method='ppr', alpha=0.15))
```

### 标志
**目的**：可扩展的 Inception 图神经网络预处理
**用例**：高效的多尺度特征
**参数**：`K`（跳数）
```python
from torch_geometric.transforms import SIGN
transform = SIGN(K=3)
```

## 特征变换

### OneHotDegree
**目的**：One-hot 编码节点度数
**用例**：度作为特征
**参数**：`max_degree`、`cat`（与现有功能连接）
```python
from torch_geometric.transforms import OneHotDegree
transform = OneHotDegree(max_degree=100)
```

### 本地学位档案
**目的**：附加本地学位简介
**用例**：结构节点特征
```python
from torch_geometric.transforms import LocalDegreeProfile
transform = LocalDegreeProfile()
```

### 常数
**目的**：为节点添加常量特征
**用例**：无特征的图表
**参数**：`value`、`cat`
```python
from torch_geometric.transforms import Constant
transform = Constant(value=1.0)
```

### 目标度数
**目的**：将入度保存为目标
**用例**：度数预测
**参数**：`norm`、`max_value`
```python
from torch_geometric.transforms import TargetIndegree
transform = TargetIndegree(norm=False)
```

### 添加RandomWalkPE
**目的**：添加随机游走位置编码
**用例**：位置信息
**参数**：`walk_length`、`attr_name`
```python
from torch_geometric.transforms import AddRandomWalkPE
transform = AddRandomWalkPE(walk_length=20)
```

### 添加拉普拉斯特征向量PE
**目的**：添加拉普拉斯特征向量位置编码
**用例**：光谱位置信息
**参数**：`k`（特征向量数量），`attr_name`
```python
from torch_geometric.transforms import AddLaplacianEigenvectorPE
transform = AddLaplacianEigenvectorPE(k=10)
```

### 添加元路径
**目的**：添加元路径引起的边缘
**用例**：异构图
**参数**：`metapaths`、`drop_orig_edges`、`drop_unconnected_nodes`
```python
from torch_geometric.transforms import AddMetaPaths
metapaths = [[('author', 'paper'), ('paper', 'author')]]  # Co-authorship
transform = AddMetaPaths(metapaths)
```

### SVDFeatureReduction
**目的**：通过 SVD 降低特征维度
**用例**：降维
**参数**：`out_channels`
```python
from torch_geometric.transforms import SVDFeatureReduction
transform = SVDFeatureReduction(out_channels=64)
```

## 视觉/空间变换

### 中心
**目的**：中心节点位置
**用例**：点云预处理
```python
from torch_geometric.transforms import Center
transform = Center()
```

### 归一化尺度
**目的**：将位置标准化为单位球体
**用例**：点云标准化
```python
from torch_geometric.transforms import NormalizeScale
transform = NormalizeScale()
```

### 标准化旋转
**目的**：旋转到主成分
**用例**：旋转不变学习
**参数**：`max_points`
```python
from torch_geometric.transforms import NormalizeRotation
transform = NormalizeRotation()
```

### 距离
**用途**：将欧氏距离保存为边缘属性
**用例**：空间图
**参数**：`norm`、`max_value`、`cat`
```python
from torch_geometric.transforms import Distance
transform = Distance(norm=False, cat=False)
```

### 笛卡尔
**用途**：将相对笛卡尔坐标保存为边缘属性
**用例**：空间关系
**参数**：`norm`、`max_value`、`cat`
```python
from torch_geometric.transforms import Cartesian
transform = Cartesian(norm=False)
```

### 极地
**用途**：将极坐标保存为边缘属性
**用例**：2D 空间图
**参数**：`norm`、`max_value`、`cat`
```python
from torch_geometric.transforms import Polar
transform = Polar(norm=False)
```

### 球形
**用途**：将球面坐标保存为边缘属性
**用例**：3D 空间图
**参数**：`norm`、`max_value`、`cat`
```python
from torch_geometric.transforms import Spherical
transform = Spherical(norm=False)
```

### 本地笛卡尔
**用途**：保存局部坐标系中的坐标
**用例**：局部空间特征
**参数**：`norm`、`cat`
```python
from torch_geometric.transforms import LocalCartesian
transform = LocalCartesian()
```

### PointPairFeatures 点对特征
**目的**：计算点对特征
**用例**：3D 配准、对应
**参数**：`cat`
```python
from torch_geometric.transforms import PointPairFeatures
transform = PointPairFeatures()
```

## 数据增强

### 随机抖动
**目的**：随机抖动节点位置
**用例**：点云增强
**参数**：`translate`、`scale`
```python
from torch_geometric.transforms import RandomJitter
transform = RandomJitter(0.01)
```

### 随机翻转
**目的**：沿轴随机翻转位置
**用例**：几何增强
**参数**：`axis`、`p`（概率）
```python
from torch_geometric.transforms import RandomFlip
transform = RandomFlip(axis=0, p=0.5)
```

### 随机尺度
**目的**：随机缩放位置
**用例**：规模扩大
**参数**：`scales`（最小值，最大值）
```python
from torch_geometric.transforms import RandomScale
transform = RandomScale((0.9, 1.1))
```

### 随机旋转
**目的**：随机轮换位置
**用例**：旋转增强
**参数**：`degrees`（范围）、`axis`（旋转轴）
```python
from torch_geometric.transforms import RandomRotate
transform = RandomRotate(degrees=15, axis=2)
```

### 随机剪切
**目的**：随机剪切位置
**用例**：几何增强
**参数**：`shear`（范围）
```python
from torch_geometric.transforms import RandomShear
transform = RandomShear(0.1)
```

### 随机翻译
**目的**：随机平移位置
**用例**：翻译增强
**参数**：`translate`（范围）
```python
from torch_geometric.transforms import RandomTranslate
transform = RandomTranslate(0.1)
```

### 线性变换
**目的**：应用线性变换矩阵
**用例**：自定义几何变换
**参数**：`matrix`
```python
from torch_geometric.transforms import LinearTransformation
import torch
matrix = torch.eye(3)
transform = LinearTransformation(matrix)
```

## 网格处理

### 样本点
**目的**：从网格中均匀采样点
**用例**：网格到点云的转换
**参数**：`num`、`remove_faces`、`include_normals`
```python
from torch_geometric.transforms import SamplePoints
transform = SamplePoints(num=1024)
```

### 生成网格法线
**用途**：生成面/顶点法线
**用例**：网格处理
```python
from torch_geometric.transforms import GenerateMeshNormals
transform = GenerateMeshNormals()
```

### 面到边
**用途**：将网格面转换为边
**用例**：网格到图形的转换
**参数**：`remove_faces`
```python
from torch_geometric.transforms import FaceToEdge
transform = FaceToEdge()
```

## 采样和分割

### 网格采样
**目的**：对体素网格中的点进行聚类
**用例**：点云下采样
**参数**：`size`（体素大小）、`start`、`end`
```python
from torch_geometric.transforms import GridSampling
transform = GridSampling(size=0.1)
```

### 定点
**目的**：采样固定数量的点
**用例**：统一点云大小
**参数**：`num`、`replace`、`allow_duplicates`
```python
from torch_geometric.transforms import FixedPoints
transform = FixedPoints(num=2048, replace=False)
```

### 随机尺度
**目的**：通过从范围内采样来随机缩放
**用例**：规模扩大（已在上面列出）

### 虚拟节点
**用途**：添加一个连接所有节点的虚拟节点
**用例**：全局信息传播
```python
from torch_geometric.transforms import VirtualNode
transform = VirtualNode()
```

## 专门的变换

### ToSLIC
**目的**：将图像转换为超像素图（SLIC算法）
**用例**：图像作为图表
**参数**：`num_segments`、`compactness`、`add_seg`、`add_img`
```python
from torch_geometric.transforms import ToSLIC
transform = ToSLIC(num_segments=75)
```

### GCNN 范数
**目的**：对边缘应用 GCN 式标准化
**用例**：GCN 预处理
**参数**：`add_self_loops`
```python
from torch_geometric.transforms import GCNNorm
transform = GCNNorm(add_self_loops=True)
```

### 拉普拉斯LambdaMax
**目的**：计算最大拉普拉斯特征值
**用例**：ChebConv 预处理
**参数**：`normalization`、`is_undirected`
```python
from torch_geometric.transforms import LaplacianLambdaMax
transform = LaplacianLambdaMax(normalization='sym')
```

### 标准化旋转
**用途**：旋转网格/点云以与主轴对齐
**用例**：规范方向
**参数**：`max_points`
```python
from torch_geometric.transforms import NormalizeRotation
transform = NormalizeRotation()
```

## 编写并应用

### 撰写
**目的**：链接多个变换
**用例**：复杂的预处理管道
```python
from torch_geometric.transforms import Compose
transform = Compose([
    Center(),
    NormalizeScale(),
    KNNGraph(k=6),
    Distance(norm=False),
])
```

### 基础变换
**用途**：自定义转换的基类
**用例**：实现自定义转换
```python
from torch_geometric.transforms import BaseTransform

class MyTransform(BaseTransform):
    def __init__(self, param):
        self.param = param

    def __call__(self, data):
        # Modify data
        data.x = data.x * self.param
        return data
```

## 常见变换组合

### 节点分类预处理
```python
transform = Compose([
    NormalizeFeatures(),
    RandomNodeSplit(num_val=0.1, num_test=0.2),
])
```

### 点云处理
```python
transform = Compose([
    Center(),
    NormalizeScale(),
    RandomRotate(degrees=15, axis=2),
    RandomJitter(0.01),
    KNNGraph(k=6),
    Distance(norm=False),
])
```

### 网格到图形
```python
transform = Compose([
    FaceToEdge(remove_faces=True),
    GenerateMeshNormals(),
    Distance(norm=True),
])
```

### 图结构增强
```python
transform = Compose([
    ToUndirected(),
    AddSelfLoops(),
    RemoveIsolatedNodes(),
    GCNNorm(),
])
```

### 异构图预处理
```python
transform = Compose([
    AddMetaPaths(metapaths=[
        [('author', 'paper'), ('paper', 'author')],
        [('author', 'paper'), ('paper', 'conference'), ('conference', 'paper'), ('paper', 'author')]
    ]),
    RandomNodeSplit(split='train_rest', num_val=0.1, num_test=0.2),
])
```

### 链接预测
```python
transform = Compose([
    NormalizeFeatures(),
    RandomLinkSplit(num_val=0.1, num_test=0.2, is_undirected=True),
])
```

## 使用提示

1. **顺序很重要**：在特征变换之前应用结构变换
2. **缓存**：某些转换（如 GDC）成本高昂——应用一次
3. **增强**：仅在训练期间使用随机*变换
4. **谨慎组合**：太多转换会减慢数据加载速度
5. **自定义转换**：继承自`BaseTransform`用于自定义逻辑
6. **预转换**：在数据集处理期间应用一次昂贵的转换：
   ```python
   dataset = MyDataset(root='/tmp', pre_transform=ExpensiveTransform())
   ```
7. **动态变换**：在训练期间应用廉价变换：
   ```python
   dataset = MyDataset(root='/tmp', transform=CheapTransform())
   ```

## 性能考虑因素

**昂贵的转换**（作为 pre_transform 应用）：
- 环球数据中心
- 标志
- KNNGraph（用于大型点云）
- 添加拉普拉斯特征向量PE
- SVD特征缩减

**便宜的变换**（作为变换应用）：
- 标准化特征
- 无向
- 添加自循环
- 随机*增强
- 到设备

**示例**：
```python
from torch_geometric.datasets import Planetoid
from torch_geometric.transforms import Compose, GDC, NormalizeFeatures

# Expensive preprocessing done once
pre_transform = GDC(
    self_loop_weight=1,
    normalization_in='sym',
    diffusion_kwargs=dict(method='ppr', alpha=0.15)
)

# Cheap transform applied each time
transform = NormalizeFeatures()

dataset = Planetoid(
    root='/tmp/Cora',
    name='Cora',
    pre_transform=pre_transform,
    transform=transform
)
```