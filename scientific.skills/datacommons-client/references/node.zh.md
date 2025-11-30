<!-- 此文件由机器翻译自 node.md -->

# 节点端点-知识图谱探索

## 目的

节点端点从数据共享知识图中检索属性关系和值。它返回有关连接节点的有向边（属性）的信息，从而能够发现图结构内的连接。

## 核心能力

Node API 执行三个主要功能：
1. 检索与节点关联的属性标签
2. 获取跨节点特定属性的值
3. 发现通过关系链接的所有连接节点

## 可用方法

### 1. fetch()

使用带有箭头表示法的关系表达式检索属性。

**关键参数：**
- `node_dcids`：目标节点标识符
- `expression`：使用箭头的关系语法（`->`、`<-`、`<-*`）
- `all_pages`：启用分页（默认值：True）
- `next_token`：继续分页结果

**箭头符号：**
- `->`：传出属性（从节点到值）
- `<-`：传入属性（从值到节点）
- `<-*`：多跳传入遍历

**用法示例：**
```python
from datacommons_client import DataCommonsClient

client = DataCommonsClient()

# Get outgoing properties from California
response = client.node.fetch(
    node_dcids=["geoId/06"],
    expression="->name"
)

# Get incoming properties (what points to this node)
response = client.node.fetch(
    node_dcids=["geoId/06"],
    expression="<-containedInPlace"
)
```

### 2. fetch_property_labels()

获取属性标签而不检索值 - 对于发现存在哪些属性很有用。

**参数：**
- `node_dcids`：节点标识符
- `out`：布尔值 - 对于传出属性为 True，对于传入属性为 False

**用法示例：**
<<<代码块_1>>>

### 3. fetch_property_values()

使用可选过滤器获取特定属性值。

**参数：**
- `node_dcids`：节点标识符
- `property`：要查询的属性名称
- `out`：方向（传出为 True，传入为 False）
- `limit`：返回值的最大数量

**用法示例：**
<<<代码块_2>>>

### 4. fetch_all_classes()

列出 Data Commons 图中的所有实体类型（类节点）。

**用法示例：**
<<<代码块_3>>>

### 5. fetch_entity_names()

按 DCID 以选定语言查找实体名称。

**参数：**
- `node_dcids`：实体标识符
- `language`：语言代码（默认值：“en”）

**用法示例：**
<<<代码块_4>>>

### 6. 放置层次结构方法

这些方法导航地理关系：

#### fetch_place_children()
获得直接儿童名额。

**用法示例：**
<<<代码块_5>>>

#### fetch_place_descendants()
检索完整的子层次结构（递归）。

**用法示例：**
<<<代码块_6>>>

#### fetch_place_parents()
获取直接父位置。

**用法示例：**
```python
# Get parent of San Francisco
parents = client.node.fetch_place_parents(
    node_dcids=["geoId/0667000"]
)
```

#### fetch_place_ancestors()
检索完整的父系谱系。

**用法示例：**
```python
# Get all ancestors of San Francisco (CA, USA, etc.)
ancestors = client.node.fetch_place_ancestors(
    node_dcids=["geoId/0667000"]
)
```

### 7. fetch_statvar_constraints()

访问统计变量的约束属性 - 对于理解变量定义和约束很有用。

**用法示例：**
```python
constraints = client.node.fetch_statvar_constraints(
    node_dcids=["Count_Person"]
)
```

## 响应格式

方法返回：
- **NodeResponse 对象** 具有 `.to_dict()`、`.to_json()` 和 `.nextToken` 属性
- **字典**用于实体名称和地点层次结构方法

## 分页

对于较大的响应：
1. 设置`all_pages=False`以块的形式接收数据
2. 响应包含 `nextToken` 值
3. 使用该令牌重新查询以获取后续页面

**示例：**
```python
# First page
response = client.node.fetch(
    node_dcids=["country/USA"],
    expression="<-containedInPlace",
    all_pages=False
)

# Get next page if available
if response.nextToken:
    next_response = client.node.fetch(
        node_dcids=["country/USA"],
        expression="<-containedInPlace",
        next_token=response.nextToken
    )
```

## 常见用例

### 用例 1：探索可用属性

```python
# Discover what properties an entity has
labels = client.node.fetch_property_labels(
    node_dcids=["geoId/06"],
    out=True
)
print(labels)  # Shows all outgoing properties like 'name', 'latitude', etc.
```

### 用例 2：导航地理层次结构

```python
# Get all counties in California
counties = client.node.fetch_place_children(
    node_dcids=["geoId/06"]
)

# Filter for specific type if needed
county_dcids = [child for child in counties["geoId/06"]
                if "County" in child]
```

### 用例 3：构建实体关系

```python
# Find all entities that reference a specific node
references = client.node.fetch(
    node_dcids=["geoId/06"],
    expression="<-location"
)
```

## 重要提示

- 首先使用`fetch_property_labels()`来发现可用的属性
- Node API 无法解析复杂的关系表达式——使用更简单的表达式或分成多个查询
- 对于具有弧关系的链接实体属性，将节点 API 查询与观察 API 结合起来
- 放置层次结构方法返回字典，而不是 NodeResponse 对象