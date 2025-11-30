<!-- 此文件由机器翻译自 api_reference.md -->

# ETE 工具包 API 参考

## 概述

ETE（树探索环境）是一个用于系统发育树操作、分析和可视化的 Python 工具包。该参考涵盖了主要的类和方法。

## 核心课程

### TreeNode（别名：树）

表示具有分层节点组织的树结构的基本类。

**构造函数：**
```python
from ete3 import Tree
t = Tree(newick=None, format=0, dist=None, support=None, name=None)
```

**参数：**
- `newick`：Newick 字符串或文件路径
- `format`：Newick 格式 (0-100)。常见格式：
  - `0`：具有分支长度和名称的灵活格式
  - `1`：具有内部节点名称
  - `2`：带有引导/支持值
  - `5`：内部节点名称和分支长度
  - `8`：所有功能（名称、距离、支持）
  - `9`：仅叶名称
  - `100`：仅拓扑
- `dist`：到父级的分支长度（默认值：1.0）
- `support`：引导/置信值（默认值：1.0）
- `name`：节点标识符

### PhyloTree

用于系统发育分析的专门类，扩展了 TreeNode。

**构造函数：**
<<<代码块_1>>>

**附加参数：**
- `alignment`：对齐文件或对齐字符串的路径
- `alg_format`：“fasta”或“phylip”
- `sp_naming_function`：从节点名称中提取物种的自定义函数

### 簇树

用于层次聚类分析的类。

**构造函数：**
<<<代码块_2>>>

**参数：**
- `text_array`：带有列标题和行名称的制表符分隔矩阵

### NCBI 分类

NCBI 分类数据库操作的类。

**构造函数：**
<<<代码块_3>>>

第一个实例化将 ~300MB NCBI 分类数据库下载到 `~/.etetoolkit/taxa.sqlite`。

## 节点属性

### 基本属性

|物业 |类型 |描述 |默认 |
|----------|------|-------------|------|
| `name` | STR |节点标识符 | 《无名》|
| `dist` |浮动|到父级的分支长度 | 1.0 |
| `support` |浮动| Bootstrap/置信值 | 1.0 |
| `up` |树节点 |父节点参考|无 |
| `children` |列表 |子节点 | []|

### 自定义功能

将任何自定义数据添加到节点：
<<<代码块_4>>>

访问功能：
<<<代码块_5>>>

## 导航和遍历

### 基本导航

<<<代码块_6>>>

### 遍历策略

```python
# Three traversal strategies
for node in tree.traverse("preorder"):    # Root → Left → Right
    print(node.name)

for node in tree.traverse("postorder"):   # Left → Right → Root
    print(node.name)

for node in tree.traverse("levelorder"):  # Level by level
    print(node.name)

# Exclude root
for node in tree.iter_descendants("postorder"):
    print(node.name)
```

### 获取节点

```python
# Get all leaves
leaves = tree.get_leaves()
for leaf in tree:  # Shortcut iteration
    print(leaf.name)

# Get all descendants
descendants = tree.get_descendants()

# Get ancestors
ancestors = node.get_ancestors()

# Get specific nodes by attribute
nodes = tree.search_nodes(name="NodeA")
node = tree & "NodeA"  # Shortcut syntax

# Get leaves by name
leaves = tree.get_leaves_by_name("LeafA")

# Get common ancestor
ancestor = tree.get_common_ancestor("LeafA", "LeafB", "LeafC")

# Custom filtering
filtered = [n for n in tree.traverse() if n.dist > 0.5 and n.is_leaf()]
```

### 迭代器方法（内存高效）

```python
# For large trees, use iterators
for match in tree.iter_search_nodes(name="X"):
    if some_condition:
        break  # Stop early

for leaf in tree.iter_leaves():
    process(leaf)

for descendant in node.iter_descendants():
    process(descendant)
```

## 树的构建和修改

### 从头开始创建树

```python
# Empty tree
t = Tree()

# Add children
child1 = t.add_child(name="A", dist=1.0)
child2 = t.add_child(name="B", dist=2.0)

# Add siblings
sister = child1.add_sister(name="C", dist=1.5)

# Populate with random topology
t.populate(10)  # Creates 10 random leaves
t.populate(5, names_library=["A", "B", "C", "D", "E"])
```

### 移除和删除节点

```python
# Detach: removes entire subtree
node.detach()
# or
parent.remove_child(node)

# Delete: removes node, reconnects children to parent
node.delete()
# or
parent.remove_child(node)
```

### 修剪

仅保留指定的叶子：
```python
# Keep only these leaves, remove all others
tree.prune(["A", "B", "C"])

# Preserve original branch lengths
tree.prune(["A", "B", "C"], preserve_branch_length=True)
```

### 树连接

```python
# Attach one tree as child of another
t1 = Tree("(A,(B,C));")
t2 = Tree("((D,E),(F,G));")
A = t1 & "A"
A.add_child(t2)
```

### 树复制

```python
# Four copy methods
copy1 = tree.copy()  # Default: cpickle (preserves types)
copy2 = tree.copy("newick")  # Fastest: basic topology
copy3 = tree.copy("newick-extended")  # Includes custom features as text
copy4 = tree.copy("deepcopy")  # Slowest: handles complex objects
```

## 树操作

### 生根

```python
# Set outgroup (reroot tree)
outgroup_node = tree & "OutgroupLeaf"
tree.set_outgroup(outgroup_node)

# Midpoint rooting
midpoint = tree.get_midpoint_outgroup()
tree.set_outgroup(midpoint)

# Unroot tree
tree.unroot()
```

### 解决多态性问题

```python
# Resolve multifurcations to bifurcations
tree.resolve_polytomy(recursive=False)  # Single node only
tree.resolve_polytomy(recursive=True)   # Entire tree
```

### 阶梯化

```python
# Sort branches by size
tree.ladderize()
tree.ladderize(direction=1)  # Ascending order
```

### 转换为超度量

```python
# Make all leaves equidistant from root
tree.convert_to_ultrametric()
tree.convert_to_ultrametric(tree_length=100)  # Specific total length
```

## 距离与比较

### 距离计算

```python
# Branch length distance between nodes
dist = tree.get_distance("A", "B")
dist = nodeA.get_distance(nodeB)

# Topology-only distance (count nodes)
dist = tree.get_distance("A", "B", topology_only=True)

# Farthest node
farthest, distance = node.get_farthest_node()
farthest_leaf, distance = node.get_farthest_leaf()
```

### 单系测试

```python
# Check if values form monophyletic group
is_mono, clade_type, base_node = tree.check_monophyly(
    values=["A", "B", "C"],
    target_attr="name"
)
# Returns: (bool, "monophyletic"|"paraphyletic"|"polyphyletic", node)

# Get all monophyletic clades
monophyletic_nodes = tree.get_monophyletic(
    values=["A", "B", "C"],
    target_attr="name"
)
```

### 树比较

```python
# Robinson-Foulds distance
rf, max_rf, common_leaves, parts_t1, parts_t2 = t1.robinson_foulds(t2)
print(f"RF distance: {rf}/{max_rf}")

# Normalized RF distance
result = t1.compare(t2)
norm_rf = result["norm_rf"]  # 0.0 to 1.0
ref_edges = result["ref_edges_in_source"]
```

## 输入/输出

### 阅读树

```python
# From string
t = Tree("(A:1,(B:1,(C:1,D:1):0.5):0.5);")

# From file
t = Tree("tree.nw")

# With format
t = Tree("tree.nw", format=1)
```

### 写树

```python
# To string
newick = tree.write()
newick = tree.write(format=1)
newick = tree.write(format=1, features=["support", "custom_feature"])

# To file
tree.write(outfile="output.nw")
tree.write(format=5, outfile="output.nw", features=["name", "dist"])

# Custom leaf function (for collapsing)
def is_leaf(node):
    return len(node) <= 3  # Treat small clades as leaves

newick = tree.write(is_leaf_fn=is_leaf)
```

### 树渲染

```python
# Show interactive GUI
tree.show()

# Render to file (PNG, PDF, SVG)
tree.render("tree.png")
tree.render("tree.pdf", w=200, units="mm")
tree.render("tree.svg", dpi=300)

# ASCII representation
print(tree)
print(tree.get_ascii(show_internal=True, compact=False))
```

## 性能优化

### 缓存内容

对于频繁访问节点内容：
```python
# Cache all node contents
node2content = tree.get_cached_content()

# Fast lookup
for node in tree.traverse():
    leaves = node2content[node]
    print(f"Node has {len(leaves)} leaves")
```

### 预计算距离

```python
# For multiple distance queries
node2dist = {}
for node in tree.traverse():
    node2dist[node] = node.get_distance(tree)
```

## PhyloTree 特定方法

### 序列比对

```python
# Link alignment
tree.link_to_alignment("alignment.fasta", alg_format="fasta")

# Access sequences
for leaf in tree:
    print(f"{leaf.name}: {leaf.sequence}")
```

### 物种命名

```python
# Default: first 3 letters
# Custom function
def get_species(node_name):
    return node_name.split("_")[0]

tree.set_species_naming_function(get_species)

# Manual setting
for leaf in tree:
    leaf.species = extract_species(leaf.name)
```

### 进化事件

```python
# Detect duplication/speciation events
events = tree.get_descendant_evol_events()

for node in tree.traverse():
    if hasattr(node, "evoltype"):
        print(f"{node.name}: {node.evoltype}")  # "D" or "S"

# With species tree
species_tree = Tree("(human, (chimp, gorilla));")
events = tree.get_descendant_evol_events(species_tree=species_tree)
```

### 基因树操作

```python
# Get species trees from duplicated gene families
species_trees = tree.get_speciation_trees()

# Split by duplication events
subtrees = tree.split_by_dups()

# Collapse lineage-specific expansions
tree.collapse_lineage_specific_expansions()
```

## NCBI 分类方法

### 数据库操作

```python
from ete3 import NCBITaxa
ncbi = NCBITaxa()

# Update database
ncbi.update_taxonomy_database()
```

### 查询分类法

```python
# Get taxid from name
taxid = ncbi.get_name_translator(["Homo sapiens"])
# Returns: {'Homo sapiens': [9606]}

# Get name from taxid
names = ncbi.get_taxid_translator([9606, 9598])
# Returns: {9606: 'Homo sapiens', 9598: 'Pan troglodytes'}

# Get rank
rank = ncbi.get_rank([9606])
# Returns: {9606: 'species'}

# Get lineage
lineage = ncbi.get_lineage(9606)
# Returns: [1, 131567, 2759, ..., 9606]

# Get descendants
descendants = ncbi.get_descendant_taxa("Primates")
descendants = ncbi.get_descendant_taxa("Primates", collapse_subspecies=True)
```

### 构建分类树

```python
# Get minimal tree connecting taxa
tree = ncbi.get_topology([9606, 9598, 9593])  # Human, chimp, gorilla

# Annotate tree with taxonomy
tree.annotate_ncbi_taxa()

# Access taxonomy info
for node in tree.traverse():
    print(f"{node.sci_name} ({node.taxid}) - Rank: {node.rank}")
```

## ClusterTree 方法

### 链接到数据

```python
# Link matrix to tree
tree.link_to_arraytable(matrix_string)

# Access profiles
for leaf in tree:
    print(leaf.profile)  # Numerical array
```

### 集群指标

```python
# Get silhouette coefficient
silhouette = tree.get_silhouette()

# Get Dunn index
dunn = tree.get_dunn()

# Inter/intra cluster distances
inter = node.intercluster_dist
intra = node.intracluster_dist

# Standard deviation
dev = node.deviation
```

### 距离指标

支持的指标：
- `"euclidean"`：欧几里德距离
- `"pearson"`：皮尔逊相关性
- `"spearman"`：斯皮尔曼等级相关性

```python
tree.dist_to(node2, metric="pearson")
```

## 常见错误处理

```python
# Check if tree is empty
if tree.children:
    print("Tree has children")

# Check if node exists
nodes = tree.search_nodes(name="X")
if nodes:
    node = nodes[0]

# Safe feature access
value = getattr(node, "feature_name", default_value)

# Check format compatibility
try:
    tree.write(format=1)
except:
    print("Tree lacks internal node names")
```

## 最佳实践

1. **使用适当的遍历**：后序为自下而上，预序为自上而下
2. **重复访问的缓存**：频繁查询时使用`get_cached_content()`
3. **对大树使用迭代器**：内存高效处理
4. **保留分支长度**：修剪时使用 `preserve_branch_length=True`
5. **明智地选择复制方法**：“newick”用于速度，“cpickle”用于完全保真度
6. **验证单系**：检查返回的进化枝类型（单系/并系/多系）
7. **使用 PhyloTree 进行系统发育**：进化分析的专门方法
8. **缓存 NCBI 查询**：存储结果以避免重复数据库访问