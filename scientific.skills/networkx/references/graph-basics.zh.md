<!-- 此文件由机器翻译自 graph-basics.md -->

# NetworkX 图基础知识

## 图表类型

NetworkX 支持四种主要图形类：

### 图（无向）
```python
import networkx as nx
G = nx.Graph()
```
- 节点之间具有单边的无向图
- 不允许有平行边
- 边缘是双向的

### 有向图（有向）
<<<代码块_1>>>
- 具有单向连接的有向图
- 边缘方向很重要：(u, v) ≠ (v, u)
- 用于建模有向关系

### MultiGraph（无向多边）
<<<代码块_2>>>
- 允许相同节点对之间存在多个边
- 对于建模多重关系很有用

### MultiDiGraph（有向多边）
<<<代码块_3>>>
- 节点之间具有多条边的有向图
- 结合了 DiGraph 和 MultiGraph 的功能

## 创建和添加节点

### 单节点添加
<<<代码块_4>>>

### 批量节点添加
<<<代码块_5>>>

### 具有属性的节点
<<<代码块_6>>>

### 重要的节点属性
- 节点可以是任何可哈希的 Python 对象：字符串、元组、数字、自定义对象
- 节点属性存储为键值对
- 为了清晰起见，使用有意义的节点标识符

## 创建和添加边

### 单边加法
```python
G.add_edge(1, 2)
G.add_edge('gene_A', 'gene_B')
```

### 批量边缘添加
```python
G.add_edges_from([(1, 2), (1, 3), (2, 4)])
G.add_edges_from(edge_list)
```

### 带属性的边
```python
G.add_edge(1, 2, weight=4.7, relation='interacts')
G.add_edges_from([
    (1, 2, {'weight': 4.7}),
    (2, 3, {'weight': 8.2, 'color': 'blue'})
])
```

### 从边缘列表添加属性
```python
# From pandas DataFrame
import pandas as pd
df = pd.DataFrame({'source': [1, 2], 'target': [2, 3], 'weight': [4.7, 8.2]})
G = nx.from_pandas_edgelist(df, 'source', 'target', edge_attr='weight')
```

## 检查图结构

### 基本属性
```python
# Get collections
G.nodes              # NodeView of all nodes
G.edges              # EdgeView of all edges
G.adj                # AdjacencyView for neighbor relationships

# Count elements
G.number_of_nodes()  # Total node count
G.number_of_edges()  # Total edge count
len(G)              # Number of nodes (shorthand)

# Degree information
G.degree()          # DegreeView of all node degrees
G.degree(1)         # Degree of specific node
list(G.degree())    # List of (node, degree) pairs
```

### 检查存在性
```python
# Check if node exists
1 in G              # Returns True/False
G.has_node(1)

# Check if edge exists
G.has_edge(1, 2)
```

### 访问邻居
```python
# Get neighbors of node 1
list(G.neighbors(1))
list(G[1])          # Dictionary-like access

# For directed graphs
list(G.predecessors(1))  # Incoming edges
list(G.successors(1))    # Outgoing edges
```

### 迭代元素
```python
# Iterate over nodes
for node in G.nodes:
    print(node, G.nodes[node])  # Access node attributes

# Iterate over edges
for u, v in G.edges:
    print(u, v, G[u][v])  # Access edge attributes

# Iterate with attributes
for node, attrs in G.nodes(data=True):
    print(node, attrs)

for u, v, attrs in G.edges(data=True):
    print(u, v, attrs)
```

## 修改图表

### 删除元素
```python
# Remove single node (also removes incident edges)
G.remove_node(1)

# Remove multiple nodes
G.remove_nodes_from([1, 2, 3])

# Remove edges
G.remove_edge(1, 2)
G.remove_edges_from([(1, 2), (2, 3)])
```

### 清除图表
```python
G.clear()           # Remove all nodes and edges
G.clear_edges()     # Remove only edges, keep nodes
```

## 属性和元数据

### 图级属性
```python
G.graph['name'] = 'Social Network'
G.graph['date'] = '2025-01-15'
print(G.graph)
```

### 节点属性
```python
# Set at creation
G.add_node(1, time='5pm', weight=0.5)

# Set after creation
G.nodes[1]['time'] = '6pm'
nx.set_node_attributes(G, {1: 'red', 2: 'blue'}, 'color')

# Get attributes
G.nodes[1]
G.nodes[1]['time']
nx.get_node_attributes(G, 'color')
```

### 边属性
```python
# Set at creation
G.add_edge(1, 2, weight=4.7, color='red')

# Set after creation
G[1][2]['weight'] = 5.0
nx.set_edge_attributes(G, {(1, 2): 10.5}, 'weight')

# Get attributes
G[1][2]
G[1][2]['weight']
G.edges[1, 2]
nx.get_edge_attributes(G, 'weight')
```

## 子图和视图

### 子图创建
```python
# Create subgraph from node list
nodes_subset = [1, 2, 3, 4]
H = G.subgraph(nodes_subset)  # Returns view (references original)

# Create independent copy
H = G.subgraph(nodes_subset).copy()

# Edge-induced subgraph
edge_subset = [(1, 2), (2, 3)]
H = G.edge_subgraph(edge_subset)
```

### 图表视图
```python
# Reverse view (for directed graphs)
G_reversed = G.reverse()

# Convert between directed/undirected
G_undirected = G.to_undirected()
G_directed = G.to_directed()
```

## 图形信息和诊断

### 基本信息
```python
print(nx.info(G))   # Summary of graph structure

# Density (ratio of actual edges to possible edges)
nx.density(G)

# Check if graph is directed
G.is_directed()

# Check if graph is multigraph
G.is_multigraph()
```

### 连接检查
```python
# For undirected graphs
nx.is_connected(G)
nx.number_connected_components(G)

# For directed graphs
nx.is_strongly_connected(G)
nx.is_weakly_connected(G)
```

## 重要考虑因素

### 浮点精度
一旦图形包含浮点数，由于精度限制，所有结果本质上都是近似的。小算术错误可能会影响算法结果，特别是在最小/最大计算中。

### 内存注意事项
每次脚本启动时，图形数据都必须加载到内存中。对于大型数据集，这可能会导致性能问题。考虑：
- 使用高效的数据格式（Python 对象的 pickle）
- 仅加载必要的子图
- 将图形数据库用于非常大的网络

### 节点和边移除行为
当删除节点时，与该节点相关的所有边也会自动删除。