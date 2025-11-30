<!-- 此文件由机器翻译自 resolve.md -->

# 解析端点 - 实体识别

## 目的

Resolve API 识别知识图中实体的数据共享 ID (DCID)。 Data Commons API 中的大多数查询都需要 DCID，因此解析通常是任何工作流程中的第一步。

## 关键能力

该端点目前支持**仅放置实体**并允许通过多种方法进行解析：
- **按名称**：使用描述性术语（例如“加利福尼亚州旧金山”）进行搜索
- **通过维基数据 ID**：使用外部标识符查找（例如，美国的“Q30”）
- **通过坐标**：通过纬度/经度定位地点
- **通过关系表达式**：使用合成属性的高级搜索

## 可用方法

### 1. fetch()

使用关系表达式进行一般解析——最灵活的方法。

**参数：**
- `nodes`：搜索词或标识符列表
- `property`：要搜索的属性（例如，“name”、“wikidataId”）

**用法示例：**
```python
from datacommons_client import DataCommonsClient

client = DataCommonsClient()

# Resolve by name
response = client.resolve.fetch(
    nodes=["California", "Texas"],
    property="name"
)
```

### 2. fetch_dcids_by_name()

具有可选类型过滤的基于名称的查找——最常用的方法。

**参数：**
- `names`：要解析的地名列表
- `entity_type`：可选类型过滤器（例如“城市”、“州”、“县”）

**返回：** `ResolveResponse` 对象，其中包含每个名称的候选对象

**用法示例：**
<<<代码块_1>>>

### 3. fetch_dcids_by_wikidata_id()

具有已知维基数据标识符的实体的维基数据 ID 解析。

**参数：**
- `wikidata_ids`：维基数据 ID 列表（例如“Q30”、“Q99”）

**用法示例：**
<<<代码块_2>>>

### 4. fetch_dcid_by_坐标()

地理坐标查找以查找特定纬度/经度坐标的位置。

**参数：**
- `latitude`：纬度坐标
- `longitude`：经度坐标

**返回：** 这些坐标处的地点的单个 DCID 字符串

**用法示例：**
<<<代码块_3>>>

## 响应结构

所有方法（`fetch_dcid_by_coordinates` 除外）都会返回一个 `ResolveResponse` 对象，其中包含：
- **节点**：提供的搜索词
- **候选**：具有可选元数据的匹配 DCID 列表
  - 每个候选可能包含 `dominantType` 字段以消除歧义
- **辅助方法**：
  - `to_dict()`：作为字典的完整响应
  - `to_json()`：JSON 字符串格式
  - `to_flat_dict()`：仅包含 DCID 的简化格式

**响应示例：**
<<<代码块_4>>>

## 常见用例

### 用例 1：在查询之前解析地名

大多数工作流程都是从将名称解析为 DCID 开始的：
<<<代码块_5>>>

### 用例 2：处理不明确的名称

当存在多个候选者时，请使用 `dominantType` 或更具体：
<<<代码块_6>>>

### 用例 3：批量解析

高效解析多个实体：
```python
places = [
    "San Francisco, CA",
    "Los Angeles, CA",
    "San Diego, CA",
    "Sacramento, CA"
]

response = client.resolve.fetch_dcids_by_name(names=places)

# Build mapping of name to DCID
name_to_dcid = {}
for name, result in response.to_dict().items():
    if result["candidates"]:
        name_to_dcid[name] = result["candidates"][0]["dcid"]
```

### 用例 4：基于坐标的查询

查找某个位置的行政地点：
```python
# User provides coordinates, find the place
latitude, longitude = 37.7749, -122.4194
dcid = client.resolve.fetch_dcid_by_coordinates(
    latitude=latitude,
    longitude=longitude
)

# Now query data for that place
response = client.observation.fetch(
    variable_dcids=["Count_Person", "MedianIncome_Household"],
    entity_dcids=[dcid],
    date="latest"
)
```

### 用例 5：外部 ID 集成

处理使用 Wikidata ID 的外部数据集时：
```python
# External dataset has Wikidata IDs
wikidata_ids = ["Q30", "Q99", "Q1384"]  # USA, California, New York

# Convert to Data Commons DCIDs
response = client.resolve.fetch_dcids_by_wikidata_id(
    wikidata_ids=wikidata_ids
)

# Extract DCIDs for further queries
dcids = []
for wid, result in response.to_dict().items():
    if result["candidates"]:
        dcids.append(result["candidates"][0]["dcid"])
```

## 重要限制

1. **仅地点实体**：Resolve API 目前仅支持地点实体（国家、州、城市、县等）。对于其他实体类型，必须通过其他方式（例如 Node API 探索）获取 DCID。

2. **无法解析链接实体属性**：对于涉及 `containedInPlace` 等关系的查询，请改用 Node API。

3. **歧义处理**：当存在多个候选者时，API 返回所有匹配项。应用程序必须根据上下文或附加过滤来决定哪个是正确的。

## 最佳实践

1. **始终首先解析名称**：切勿采用 DCID 格式 - 始终使用 Resolve API
2. **缓存解析**：如果重复查询同一个地方，缓存名称→DCID映射
3. **处理歧义**：检查多个候选者并使用 `entity_type` 过滤或更具体的名称
4. **验证结果**：在访问 DCID 之前，始终检查 `candidates` 列表不为空
5. **使用适当的方法**：
   - 名称 → `fetch_dcids_by_name()`
   - 坐标 → `fetch_dcid_by_coordinates()`
   - 维基数据 ID → `fetch_dcids_by_wikidata_id()`