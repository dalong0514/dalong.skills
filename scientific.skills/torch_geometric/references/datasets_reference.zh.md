<!-- 此文件由机器翻译自 datasets_reference.md -->

# PyTorch 几何数据集参考

本文档提供了 `torch_geometric.datasets` 中可用的所有数据集的综合目录。

## 引文网络

### 小行星
**用途**：节点分类、半监督学习
**网络**：Cora、CiteSeer、PubMed
**描述**：引文网络，其中节点是论文，边是引文
- **Cora**：2,708 个节点，5,429 个边，7 个类，1,433 个特征
- **CiteSeer**：3,327 个节点，4,732 个边，6 个类，3,703 个特征
- **PubMed**：19,717 个节点，44,338 个边，3 个类，500 个特征

```python
from torch_geometric.datasets import Planetoid
dataset = Planetoid(root='/tmp/Cora', name='Cora')
```

### 共同作者
**用途**：协作网络上的节点分类
**网络**：计算机科学、物理
**描述**：来自 Microsoft Academy Graph 的合着网络
- **CS**：18,333 个节点，81,894 个边，15 个类（计算机科学）
- **物理**：34,493 个节点，247,962 个边，5 个类（物理）

<<<代码块_1>>>

### 亚马逊
**用途**：产品网络上的节点分类
**网络**：计算机、照片
**描述**：亚马逊共购网络，其中节点是产品
- **计算机**：13,752 个节点，245,861 个边，10 个类
- **照片**：7,650 个节点，119,081 条边，8 个类

<<<代码块_2>>>

### 引文完整
**用途**：引文网络分析
**网络**：Cora、Cora_ML、DBLP、PubMed
**描述**：无需抽样的完整引用网络

<<<代码块_3>>>

## 图分类

### TU数据集
**用途**：图分类、图内核基准测试
**描述**：120+图分类数据集的集合
- **MUTAG**：188 个图表，2 个类别（分子化合物）
- **蛋白质**：1,113 个图表，2 个类别（蛋白质结构）
- **酶**：600 个图表，6 个类别（蛋白质酶）
- **IMDB-BINARY**：1,000 个图表，2 个类别（社交网络）
- **REDDIT-BINARY**：2,000 个图表，2 个课程（讨论线程）
- **COLLAB**：5,000 个图表，3 个类别（科学合作）
- **NCI1**：4,110 个图表，2 个类别（化合物）
- **DD**：1,178 个图表，2 个类别（蛋白质结构）

<<<代码块_4>>>

### 分子网
**用途**：分子性质预测
**数据集**：超过 10 个分子基准数据集
**描述**：全面的分子机器学习基准
- **ESOL**：水溶性（回归）
- **FreeSolv**：水合自由能（回归）
- **亲脂性**：辛醇/水分布（回归）
- **BACE**：结合结果（分类）
- **BBBP**：血脑屏障穿透（分类）
- **HIV**：HIV抑制（分类）
- **Tox21**：毒性预测（多任务分类）
- **ToxCast**：毒理学预测（多任务分类）
- **SIDER**：副作用（多任务分类）
- **ClinTox**：临床试验毒性（多任务分类）

<<<代码块_5>>>

## 分子和化学数据集

### QM7b
**用途**：分子性质预测（量子力学）
**描述**：7,211 个分子，最多有 7 个重原子
- 特性：雾化能、电子特性

<<<代码块_6>>>

### QM9
**用途**：分子性质预测（量子力学）
**描述**：约 130,000 个分子，最多含 9 个重原子（C、O、N、F）
- 性质：19种量子化学性质，包括HOMO、LUMO、能隙、能量

```python
from torch_geometric.datasets import QM9
dataset = QM9(root='/tmp/QM9')
```

### 锌
**用途**：分子生成、性质预测
**描述**：约 250,000 个类药物分子图
- 性质：受限溶解度、分子量

```python
from torch_geometric.datasets import ZINC
dataset = ZINC(root='/tmp/ZINC', subset=True)
```

### AQSOL
**用途**：水溶性预测
**描述**：约 10,000 个分子的溶解度测量

```python
from torch_geometric.datasets import AQSOL
dataset = AQSOL(root='/tmp/AQSOL')
```

### MD17
**用途**：分子动力学、力场学习
**描述**：小分子的分子动力学轨迹
- 分子：苯、尿嘧啶、萘、阿司匹林、水杨酸等。

```python
from torch_geometric.datasets import MD17
dataset = MD17(root='/tmp/MD17', name='benzene')
```

### PCQM4Mv2
**用途**：大规模分子特性预测
**描述**：来自 PubChem 的 3.8M 分子用于量子化学
- OGB 大规模挑战赛的一部分

```python
from torch_geometric.datasets import PCQM4Mv2
dataset = PCQM4Mv2(root='/tmp/PCQM4Mv2')
```

## 社交网络

### 红迪网
**用途**：大规模节点分类
**描述**：2014 年 9 月的 Reddit 帖子
- 232,965 个节点，11,606,919 个边，41 个类
- 特点：帖子内容的 TF-IDF

```python
from torch_geometric.datasets import Reddit
dataset = Reddit(root='/tmp/Reddit')
```

### Reddit2
**用途**：大规模节点分类
**描述**：更新了 Reddit 数据集并添加了更多帖子

```python
from torch_geometric.datasets import Reddit2
dataset = Reddit2(root='/tmp/Reddit2')
```

### 抽搐
**用途**：节点分类、社交网络分析
**网络**：德语、英语、西班牙语、法语、葡萄牙语、俄语
**描述**：按语言划分 Twitch 用户网络

```python
from torch_geometric.datasets import Twitch
dataset = Twitch(root='/tmp/Twitch', name='DE')
```

### 脸书
**用途**：社交网络分析、节点分类
**描述**：Facebook 页面-页面网络

```python
from torch_geometric.datasets import FacebookPagePage
dataset = FacebookPagePage(root='/tmp/Facebook')
```

### GitHub
**用途**：社交网络分析
**描述**：GitHub 开发者网络

```python
from torch_geometric.datasets import GitHub
dataset = GitHub(root='/tmp/GitHub')
```

## 知识图谱

### 实体
**用途**：链接预测、知识图嵌入
**数据集**：AIFB、MUTAG、BGS、AM
**描述**：具有类型化关系的RDF知识图

```python
from torch_geometric.datasets import Entities
dataset = Entities(root='/tmp/AIFB', name='AIFB')
```

### WordNet18
**用途**：语义网络上的链接预测
**描述**：具有 18 个关系的 WordNet 子集
- 40,943 个实体，151,442 个三元组

```python
from torch_geometric.datasets import WordNet18
dataset = WordNet18(root='/tmp/WordNet18')
```

### WordNet18RR
**用途**：链接预测（无逆关系）
**描述**：没有逆关系的精炼版本

```python
from torch_geometric.datasets import WordNet18RR
dataset = WordNet18RR(root='/tmp/WordNet18RR')
```

### FB15k-237
**用法**：Freebase 上的链接预测
**描述**：具有 237 个关系的 Freebase 子集
- 14,541 个实体，310,116 个三元组

```python
from torch_geometric.datasets import FB15k_237
dataset = FB15k_237(root='/tmp/FB15k')
```

## 异构图

### OGB_MAG
**用途**：异构图学习、节点分类
**描述**：具有多种节点/边类型的 Microsoft Academy Graph
- 节点类型：论文、作者、机构、研究领域
- 1M+ 节点，21M+ 边

```python
from torch_geometric.datasets import OGB_MAG
dataset = OGB_MAG(root='/tmp/OGB_MAG')
```

### 电影镜头
**用途**：推荐系统、链接预测
**版本**：100K、1M、10M、20M
**描述**：用户电影评级网络
- 节点类型：用户、电影
- 边缘类型：费率

```python
from torch_geometric.datasets import MovieLens
dataset = MovieLens(root='/tmp/MovieLens', model_name='100k')
```

### 互联网数据库
**用途**：异构图学习
**描述**：IMDB 电影网络
- 节点类型：电影、演员、导演

```python
from torch_geometric.datasets import IMDB
dataset = IMDB(root='/tmp/IMDB')
```

### DBLP
**用途**：异构图学习、节点分类
**描述**：DBLP参考书目网络
- 节点类型：作者、论文、术语、会议

```python
from torch_geometric.datasets import DBLP
dataset = DBLP(root='/tmp/DBLP')
```

### 最后FM
**用途**：异构推荐
**描述**：LastFM 音乐网络
- 节点类型：用户、艺术家、标签

```python
from torch_geometric.datasets import LastFM
dataset = LastFM(root='/tmp/LastFM')
```

## 时间图

### 比特币场外交易
**用途**：时间链接预测、信任网络
**描述**：随着时间的推移，比特币场外交易信任网络

```python
from torch_geometric.datasets import BitcoinOTC
dataset = BitcoinOTC(root='/tmp/BitcoinOTC')
```

### ICEWS18
**用法**：时态知识图补全
**描述**：综合危机预警系统事件

```python
from torch_geometric.datasets import ICEWS18
dataset = ICEWS18(root='/tmp/ICEWS18')
```

### GDELT
**用途**：时间事件预测
**描述**：事件、语言和语气的全球数据库

```python
from torch_geometric.datasets import GDELT
dataset = GDELT(root='/tmp/GDELT')
```

### JODIE数据集
**用途**：动态图学习
**数据集**：Reddit、维基百科、MOOC、LastFM
**描述**：时间交互网络

```python
from torch_geometric.datasets import JODIEDataset
dataset = JODIEDataset(root='/tmp/JODIE', name='Reddit')
```

## 3D 网格和点云

### 形状网
**用途**：3D 形状分类和分割
**描述**：大型3D CAD模型数据集
- 16 个类别的 16,881 个模型
- 部件级分割标签

```python
from torch_geometric.datasets import ShapeNet
dataset = ShapeNet(root='/tmp/ShapeNet', categories=['Airplane'])
```

### 模型网
**用途**：3D形状分类
**版本**：ModelNet10、ModelNet40
**描述**：用于 3D 对象分类的 CAD 模型
- ModelNet10：4,899 个模型，10 个类别
- ModelNet40：12,311 个模型，40 个类别

```python
from torch_geometric.datasets import ModelNet
dataset = ModelNet(root='/tmp/ModelNet', name='10')
```

### 浮士德
**用途**：3D形状匹配、对应
**描述**：人体扫描进行形状分析
- 10 个人 10 个姿势的 100 个网格

```python
from torch_geometric.datasets import FAUST
dataset = FAUST(root='/tmp/FAUST')
```

### 科玛
**用途**：3D网格变形
**描述**：面部表情网格
- 20,466 个带表情的 3D 面部扫描

```python
from torch_geometric.datasets import CoMA
dataset = CoMA(root='/tmp/CoMA')
```

### S3DIS
**用途**：3D语义分割
**描述**：斯坦福大学大型3D室内空间
- 6个区域，271个房间，点云数据

```python
from torch_geometric.datasets import S3DIS
dataset = S3DIS(root='/tmp/S3DIS', test_area=6)
```

## 图像和视觉数据集

### MNISTSuperpixels
**用途**：基于图的图像分类
**描述**：MNIST 图像作为超像素图
- 70,000 个图表（60k 训练，10k 测试）

```python
from torch_geometric.datasets import MNISTSuperpixels
dataset = MNISTSuperpixels(root='/tmp/MNIST')
```

### Flickr
**用途**：图像描述、节点分类
**描述**：Flickr 图片网络
- 89,250 个节点，899,756 个边

```python
from torch_geometric.datasets import Flickr
dataset = Flickr(root='/tmp/Flickr')
```

### 生产者价格指数
**用途**：蛋白质-蛋白质相互作用预测
**描述**：多图蛋白质相互作用网络
- 24 个图表，总共 2,373 个节点

```python
from torch_geometric.datasets import PPI
dataset = PPI(root='/tmp/PPI', split='train')
```

## 小经典图

### 空手道俱乐部
**用途**：社区检测、可视化
**描述**：Zachary 的空手道俱乐部网络
- 34 个节点，78 个边，2 个社区

```python
from torch_geometric.datasets import KarateClub
dataset = KarateClub()
```

## 开放图基准（OGB）

PyG 与 OGB 数据集无缝集成：

### 节点属性预测
- **ogbn-products**：亚马逊产品网络（240 万个节点）
- **ogbn-蛋白质**：蛋白质关联网络（132K 节点）
- **ogbn-arxiv**：引文网络（169K 节点）
- **ogbn-papers100M**：大型引用网络（111M 节点）
- **ogbn-mag**：异构学术图谱

### 链接属性预测
- **ogbl-ppa**：蛋白质关联网络
- **ogbl-collab**：协作网络
- **ogbl-ddi**：药物-药物相互作用网络
- **ogbl-itation2**：引文网络
- **ogbl-wikikg2**：维基数据知识图

### 图属性预测
- **ogbg-molhiv**：HIV 分子活性预测
- **ogbg-molpcba**：分子生物测定（多任务）
- **ogbg-ppa**：蛋白质功能预测
- **ogbg-code2**：代码抽象语法树

```python
from torch_geometric.datasets import OGB_MAG, OGB_PPA
# or
from ogb.nodeproppred import PygNodePropPredDataset
dataset = PygNodePropPredDataset(name='ogbn-arxiv')
```

## 综合数据集

### 假数据集
**用途**：测试、调试
**描述**：生成随机图形数据

```python
from torch_geometric.datasets import FakeDataset
dataset = FakeDataset(num_graphs=100, avg_num_nodes=50)
```

### 随机块模型数据集
**用途**：社区检测基准
**描述**：从随机块模型生成的图

```python
from torch_geometric.datasets import StochasticBlockModelDataset
dataset = StochasticBlockModelDataset(root='/tmp/SBM', num_graphs=1000)
```

### 解释器数据集
**用法**：测试可解释性方法
**描述**：具有已知解释基本事实的合成图

```python
from torch_geometric.datasets import ExplainerDataset
dataset = ExplainerDataset(num_graphs=1000)
```

## 材料科学

### QM8
**用途**：分子性质预测
**描述**：小分子的电子特性

```python
from torch_geometric.datasets import QM8
dataset = QM8(root='/tmp/QM8')
```

## 生物网络

### PPI（蛋白质-蛋白质相互作用）
已在上面的图像和视觉数据集下列出

### 字符串
**用途**：蛋白质相互作用网络
**描述**：已知和预测的蛋白质-蛋白质相互作用

```python
# Available through external sources or custom loading
```

## 使用提示

1. **从小数据集开始**：使用 Cora、KarateClub 或 ENZYMES 进行原型设计
2. **引文网络**：Planetoid数据集非常适合节点分类
3. **图分类**：TUDataset提供多样化的基准
4. **分子**：用于化学应用的 QM9、ZINC、MoleculeNet
5. **大规模**：通过 NeighborLoader 使用 Reddit、OGB 数据集
6. **异构**：OGB_MAG、MovieLens、IMDB 用于多类型图
7. **Temporal**：用于动态图学习的 JODIE、ICEWS
8. **3D**：用于几何学习的ShapeNet、ModelNet、S3DIS

## 常见模式

### 使用变换加载
```python
from torch_geometric.datasets import Planetoid
from torch_geometric.transforms import NormalizeFeatures

dataset = Planetoid(root='/tmp/Cora', name='Cora',
                    transform=NormalizeFeatures())
```

### 训练/验证/测试分割
```python
# For datasets with pre-defined splits
data = dataset[0]
train_data = data[data.train_mask]
val_data = data[data.val_mask]
test_data = data[data.test_mask]

# For graph classification
from torch_geometric.loader import DataLoader
train_dataset = dataset[:int(len(dataset) * 0.8)]
test_dataset = dataset[int(len(dataset) * 0.8):]
train_loader = DataLoader(train_dataset, batch_size=32)
```

### 自定义数据加载
```python
from torch_geometric.data import Data, Dataset

class MyCustomDataset(Dataset):
    def __init__(self, root, transform=None):
        super().__init__(root, transform)
        # Your initialization

    def len(self):
        return len(self.data_list)

    def get(self, idx):
        # Load and return data object
        return self.data_list[idx]
```