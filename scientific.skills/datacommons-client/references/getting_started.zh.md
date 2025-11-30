<!-- 此文件由机器翻译自 getting_started.md -->

# 数据共享入门

## 快速入门指南

本指南提供常见 Data Commons 工作流程的端到端示例。

## 安装和设置

```bash
# Install with Pandas support
pip install "datacommons-client[Pandas]"

# Set up API key for datacommons.org
export DC_API_KEY="your_api_key_here"
```

请求 API 密钥：https://apikeys.datacommons.org/

## 示例1：基本人口查询

查询特定地点的当前人口数量：

<<<代码块_1>>>

## 示例 2：时间序列分析

检索并绘制历史失业率：

<<<代码块_2>>>

## 示例3：地理层次结构查询

获取一个州内所有县的数据：

<<<代码块_3>>>

## 示例 4：多变量比较

比较各个实体的多个统计数据：

<<<代码块_4>>>

## 示例 5：基于坐标的查询

通过坐标查找和查询某个位置的数据：

<<<代码块_5>>>

## 示例6：数据源过滤

从特定来源查询数据的一致性：

<<<代码块_6>>>

## 示例 7：探索知识图

发现实体属性和关系：

```python
client = DataCommonsClient()

# Step 1: Explore what properties exist
entity = "geoId/06"  # California

# Get outgoing properties
out_props = client.node.fetch_property_labels(
    node_dcids=[entity],
    out=True
)

print(f"Outgoing properties for California:")
print(out_props[entity])

# Get incoming properties
in_props = client.node.fetch_property_labels(
    node_dcids=[entity],
    out=False
)

print(f"\nIncoming properties for California:")
print(in_props[entity])

# Step 2: Get specific property values
name_response = client.node.fetch_property_values(
    node_dcids=[entity],
    property="name",
    out=True
)

print(f"\nName property value:")
print(name_response.to_dict())

# Step 3: Explore hierarchy
children = client.node.fetch_place_children(node_dcids=[entity])
print(f"\nNumber of child places: {len(children[entity])}")

# Get names for first 5 children
if children[entity]:
    child_sample = children[entity][:5]
    child_names = client.node.fetch_entity_names(node_dcids=child_sample)
    print("\nSample child places:")
    for dcid, name in child_names.items():
        print(f"  {name}")
```

## 示例 8：批处理多个查询

高效查询多个实体的数据：

```python
import pandas as pd

client = DataCommonsClient()

# List of cities to analyze
cities = [
    "San Francisco, CA",
    "Los Angeles, CA",
    "San Diego, CA",
    "Sacramento, CA",
    "San Jose, CA"
]

# Resolve all cities
resolve_response = client.resolve.fetch_dcids_by_name(
    names=cities,
    entity_type="City"
)

# Build mapping
city_dcids = []
dcid_to_name = {}
for name, result in resolve_response.to_dict().items():
    if result["candidates"]:
        dcid = result["candidates"][0]["dcid"]
        city_dcids.append(dcid)
        dcid_to_name[dcid] = name

# Query multiple variables at once
variables = [
    "Count_Person",
    "Median_Income_Household",
    "UnemploymentRate_Person"
]

response = client.observation.fetch(
    variable_dcids=variables,
    entity_dcids=city_dcids,
    date="latest"
)

# Process into a comparison table
df = response.to_observations_as_records()
df['city'] = df['entity'].map(dcid_to_name)

# Create comparison table
comparison = df.pivot_table(
    values='value',
    index='city',
    columns='variable',
    aggfunc='first'
)

print("\nCalifornia Cities Comparison:")
print(comparison.to_string())

# Export to CSV
comparison.to_csv('ca_cities_comparison.csv')
print("\nData exported to ca_cities_comparison.csv")
```

## 常见模式总结

### 模式 1：名称 → DCID → 数据
```python
names = ["California"]
dcids = resolve_names(names)
data = query_observations(dcids, variables)
```

### 模式 2：坐标 → DCID → 数据
```python
dcid = resolve_coordinates(lat, lon)
data = query_observations([dcid], variables)
```

### 模式 3：父项 → 子项 → 数据
```python
children = get_place_children(parent_dcid)
data = query_observations(children, variables)
```

### 模式 4：探索→选择→查询
```python
available_vars = check_available_variables(dcids)
selected_vars = filter_relevant(available_vars)
data = query_observations(dcids, selected_vars)
```

## 错误处理最佳实践

```python
client = DataCommonsClient()

# Always check for candidates
resolve_response = client.resolve.fetch_dcids_by_name(names=["Unknown Place"])
result = resolve_response.to_dict()["Unknown Place"]

if not result["candidates"]:
    print("No matches found - try a more specific name")
    # Handle error appropriately
else:
    dcid = result["candidates"][0]["dcid"]
    # Proceed with query

# Check for multiple candidates (ambiguity)
if len(result["candidates"]) > 1:
    print(f"Multiple matches found: {len(result['candidates'])}")
    for candidate in result["candidates"]:
        print(f"  {candidate['dcid']} ({candidate.get('dominantType', 'N/A')})")
    # Let user select or use additional filtering
```

## 后续步骤

1. 探索可用的统计变量：https://datacommons.org/tools/statvar
2.浏览知识图谱：https://datacommons.org/browser/
3.阅读`references/`目录中的详细端点文档
4.查看官方文档：https://docs.datacommons.org/api/python/v2/