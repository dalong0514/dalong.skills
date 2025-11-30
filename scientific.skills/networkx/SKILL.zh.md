<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：网络x
描述：用于在 Python 中创建、分析和可视化复杂网络和图形的综合工具包。在处理网络/图形数据结构、分析实体之间的关系、计算图形算法（最短路径、中心性、聚类）、检测社区、生成合成网络或可视化网络拓扑时使用。适用于社交网络、生物网络、交通系统、引文网络以及任何涉及成对关系的领域。
---

# 网络X

## 概述

NetworkX 是一个用于创建、操作和分析复杂网络和图形的 Python 包。在处理网络或图形数据结构时使用此技能，包括社交网络、生物网络、运输系统、引文网络、知识图或涉及实体之间关系的任何系统。

## 何时使用此技能

当任务涉及以下内容时调用此技能：

- **创建图表**：从数据构建网络结构，添加具有属性的节点和边
- **图分析**：计算中心性度量、寻找最短路径、检测社区、测量聚类
- **图算法**：运行 Dijkstra、PageRank、最小生成树、最大流量等标准算法
- **网络生成**：创建用于测试或模拟的合成网络（随机、无标度、小世界模型）
- **图形 I/O**：读取或写入各种格式（边列表、GraphML、JSON、CSV、邻接矩阵）
- **可视化**：使用 matplotlib 或交互式库绘制和自定义网络可视化
- **网络比较**：检查同构、计算图指标、分析结构属性

## 核心能力

### 1. 图形创建和操作

NetworkX 支持四种主要图形类型：
- **图**：具有单边的无向图
- **有向图**：具有单向连接的有向图
- **MultiGraph**：无向图允许节点之间有多个边
- **MultiDiGraph**：具有多条边的有向图

通过以下方式创建图表：
```python
import networkx as nx

# Create empty graph
G = nx.Graph()

# Add nodes (can be any hashable type)
G.add_node(1)
G.add_nodes_from([2, 3, 4])
G.add_node("protein_A", type='enzyme', weight=1.5)

# Add edges
G.add_edge(1, 2)
G.add_edges_from([(1, 3), (2, 4)])
G.add_edge(1, 4, weight=0.8, relation='interacts')
```

**参考**：有关创建、修改、检查和管理图结构（包括使用属性和子图）的综合指南，请参阅 `references/graph-basics.md`。

### 2.图算法

NetworkX 提供了广泛的网络分析算法：

**最短路径**：
<<<代码块_1>>>

**中心性措施**：
<<<代码块_2>>>

**社区检测**：
<<<代码块_3>>>

**连接性**：
<<<代码块_4>>>

**参考**：请参阅 `references/algorithms.md` 有关所有可用算法的详细文档，包括最短路径、中心性度量、聚类、社区检测、流、匹配、树算法和图遍历。

### 3.图形生成器

创建用于测试、模拟或建模的合成网络：

**经典图表**：
<<<代码块_5>>>

**随机网络**：
<<<代码块_6>>>

**结构化网络**：
```python
# Grid graph
G = nx.grid_2d_graph(m=5, n=7)

# Random tree
G = nx.random_tree(n=100, seed=42)
```

**参考**：请参阅`references/generators.md`，了解所有图形生成器的全面覆盖，包括经典、随机、格子、二分和具有详细参数和用例的专用网络模型。

### 4. 读写图表

NetworkX 支持多种文件格式和数据源：

**文件格式**：
```python
# Edge list
G = nx.read_edgelist('graph.edgelist')
nx.write_edgelist(G, 'graph.edgelist')

# GraphML (preserves attributes)
G = nx.read_graphml('graph.graphml')
nx.write_graphml(G, 'graph.graphml')

# GML
G = nx.read_gml('graph.gml')
nx.write_gml(G, 'graph.gml')

# JSON
data = nx.node_link_data(G)
G = nx.node_link_graph(data)
```

**熊猫集成**：
```python
import pandas as pd

# From DataFrame
df = pd.DataFrame({'source': [1, 2, 3], 'target': [2, 3, 4], 'weight': [0.5, 1.0, 0.75]})
G = nx.from_pandas_edgelist(df, 'source', 'target', edge_attr='weight')

# To DataFrame
df = nx.to_pandas_edgelist(G)
```

**矩阵格式**：
```python
import numpy as np

# Adjacency matrix
A = nx.to_numpy_array(G)
G = nx.from_numpy_array(A)

# Sparse matrix
A = nx.to_scipy_sparse_array(G)
G = nx.from_scipy_sparse_array(A)
```

**参考**：有关所有 I/O 格式（包括 CSV、SQL 数据库、Cytoscape、DOT）的完整文档，以及不同用例的格式选择指南，请参阅 `references/io.md`。

### 5.可视化

创建清晰且信息丰富的网络可视化：

**基本可视化**：
```python
import matplotlib.pyplot as plt

# Simple draw
nx.draw(G, with_labels=True)
plt.show()

# With layout
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos=pos, with_labels=True, node_color='lightblue', node_size=500)
plt.show()
```

**定制**：
```python
# Color by degree
node_colors = [G.degree(n) for n in G.nodes()]
nx.draw(G, node_color=node_colors, cmap=plt.cm.viridis)

# Size by centrality
centrality = nx.betweenness_centrality(G)
node_sizes = [3000 * centrality[n] for n in G.nodes()]
nx.draw(G, node_size=node_sizes)

# Edge weights
edge_widths = [3 * G[u][v].get('weight', 1) for u, v in G.edges()]
nx.draw(G, width=edge_widths)
```

**布局算法**：
```python
# Spring layout (force-directed)
pos = nx.spring_layout(G, seed=42)

# Circular layout
pos = nx.circular_layout(G)

# Kamada-Kawai layout
pos = nx.kamada_kawai_layout(G)

# Spectral layout
pos = nx.spectral_layout(G)
```

**出版质量**：
```python
plt.figure(figsize=(12, 8))
pos = nx.spring_layout(G, seed=42)
nx.draw(G, pos=pos, node_color='lightblue', node_size=500,
        edge_color='gray', with_labels=True, font_size=10)
plt.title('Network Visualization', fontsize=16)
plt.axis('off')
plt.tight_layout()
plt.savefig('network.png', dpi=300, bbox_inches='tight')
plt.savefig('network.pdf', bbox_inches='tight')  # Vector format
```

**参考**：请参阅 `references/visualization.md` 有关可视化技术的大量文档，包括布局算法、自定义选项、使用 Plotly 和 PyVis 进行交互式可视化、3D 网络以及出版物质量的图形创建。

## 使用 NetworkX

### 安装

确保 NetworkX 已安装：
```python
# Check if installed
import networkx as nx
print(nx.__version__)

# Install if needed (via bash)
# uv pip install networkx
# uv pip install networkx[default]  # With optional dependencies
```

### 常见工作流程模式

大多数 NetworkX 任务都遵循以下模式：

1. **创建或加载图表**：
   ```python
   # From scratch
   G = nx.Graph()
   G.add_edges_from([(1, 2), (2, 3), (3, 4)])

   # Or load from file/data
   G = nx.read_edgelist('data.txt')
   ```

2. **检查结构**：
   ```python
   print(f"Nodes: {G.number_of_nodes()}")
   print(f"Edges: {G.number_of_edges()}")
   print(f"Density: {nx.density(G)}")
   print(f"Connected: {nx.is_connected(G)}")
   ```

3. **分析**：
   ```python
   # Compute metrics
   degree_cent = nx.degree_centrality(G)
   avg_clustering = nx.average_clustering(G)

   # Find paths
   path = nx.shortest_path(G, source=1, target=4)

   # Detect communities
   communities = community.greedy_modularity_communities(G)
   ```

4. **可视化**：
   ```python
   pos = nx.spring_layout(G, seed=42)
   nx.draw(G, pos=pos, with_labels=True)
   plt.show()
   ```
5. **导出结果**：
   ```python
   # Save graph
   nx.write_graphml(G, 'analyzed_network.graphml')

   # Save metrics
   df = pd.DataFrame({
       'node': list(degree_cent.keys()),
       'centrality': list(degree_cent.values())
   })
   df.to_csv('centrality_results.csv', index=False)
   ```

### 重要考虑因素

**浮点精度**：当图形包含浮点数时，由于精度限制，所有结果本质上都是近似值。这可能会影响算法结果，特别是在最小/最大计算中。

**内存和性能**：每次运行脚本时，图形数据都必须加载到内存中。对于大型网络：
- 使用适当的数据结构（大型稀疏图的稀疏矩阵）
- 考虑仅加载必要的子图
- 使用高效的文件格式（Python 对象的 pickle、压缩格式）
- 利用非常大的网络的近似算法（例如，中心性计算中的 `k` 参数）

**节点和边类型**：
- 节点可以是任何可哈希的 Python 对象（数字、字符串、元组、自定义对象）
- 为了清晰起见，使用有意义的标识符
- 删除节点时，所有关联边都会自动删除

**随机种子**：始终设置随机种子，以实现随机图生成和力导向布局中的再现性：
```python
G = nx.erdos_renyi_graph(n=100, p=0.1, seed=42)
pos = nx.spring_layout(G, seed=42)
```

## 快速参考

### 基本操作
```python
# Create
G = nx.Graph()
G.add_edge(1, 2)

# Query
G.number_of_nodes()
G.number_of_edges()
G.degree(1)
list(G.neighbors(1))

# Check
G.has_node(1)
G.has_edge(1, 2)
nx.is_connected(G)

# Modify
G.remove_node(1)
G.remove_edge(1, 2)
G.clear()
```

### 基本算法
```python
# Paths
nx.shortest_path(G, source, target)
nx.all_pairs_shortest_path(G)

# Centrality
nx.degree_centrality(G)
nx.betweenness_centrality(G)
nx.closeness_centrality(G)
nx.pagerank(G)

# Clustering
nx.clustering(G)
nx.average_clustering(G)

# Components
nx.connected_components(G)
nx.strongly_connected_components(G)  # Directed

# Community
community.greedy_modularity_communities(G)
```

### 文件 I/O 快速参考
```python
# Read
nx.read_edgelist('file.txt')
nx.read_graphml('file.graphml')
nx.read_gml('file.gml')

# Write
nx.write_edgelist(G, 'file.txt')
nx.write_graphml(G, 'file.graphml')
nx.write_gml(G, 'file.gml')

# Pandas
nx.from_pandas_edgelist(df, 'source', 'target')
nx.to_pandas_edgelist(G)
```

## 资源

该技能包括全面的参考文档：

### 参考文献/graph-basics.md
有关图类型、创建和修改图、添加节点和边、管理属性、检查结构以及使用子图的详细指南。

### 参考文献/algorithms.md
完整覆盖 NetworkX 算法，包括最短路径、中心性度量、连通性、聚类、社区检测、流算法、树算法、匹配、着色、同构和图遍历。

### 参考文献/generators.md
有关图生成器的综合文档，包括经典图、随机模型（Erdős-Rényi、Barabási-Albert、Watts-Strogatz）、格子、树、社交网络模型和专用生成器。

### 参考文献/io.md
读写各种格式图形的完整指南：边列表、邻接列表、GraphML、GML、JSON、CSV、Pandas DataFrames、NumPy 数组、SciPy 稀疏矩阵、数据库集成和格式选择指南。

### 参考文献/可视化.md
有关可视化技术的大量文档，包括布局算法、自定义节点和边缘外观、标签、使用 Plotly 和 PyVis 进行交互式可视化、3D 网络、二分布局以及创建出版物质量的图形。

## 其他资源

- **官方文档**：https://networkx.org/documentation/latest/
- **教程**：https://networkx.org/documentation/latest/tutorial.html
- **图库**：https://networkx.org/documentation/latest/auto_examples/index.html
- **GitHub**：https://github.com/networkx/networkx