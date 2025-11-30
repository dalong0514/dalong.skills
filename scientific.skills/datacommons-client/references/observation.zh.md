<!-- 此文件由机器翻译自 observation.md -->

# 观察端点-统计数据查询

## 目的

观察 API 检索统计观察结果——链接实体、变量和特定日期的数据点。示例包括：
- 《2020 年美国人口》
- “加州随时间变化的 GDP”
- “一个州所有县的失业率”

## 核心方法

### 1. fetch()

使用灵活的实体规范检索观测值的主要方法。

**关键参数：**
- `variable_dcids`（必需）：统计变量标识符列表
- `entity_dcids` 或 `entity_expression`（必需）：通过 ID 或关系表达式指定实体
- `date`（可选）：默认为“最新”。接受：
  - ISO-8601 格式（例如“2020”、“2020-01”、“2020-01-15”）
  - “全部”表示完整的时间序列
  - “最新”表示最新数据
- `select`（可选）：控制返回字段
  - 默认：`["date", "entity", "variable", "value"]`
  - 替代方案：`["entity", "variable", "facet"]` 在没有数据的情况下检查可用性
- `filter_facet_domains`：按数据源域过滤
- `filter_facet_ids`：按特定方面 ID 过滤

**响应结构：**
数据按变量→实体分层组织，有关“方面”（数据源）的元数据包括：
- 出处 URL
- 测量方法
- 观察期
- 导入名称

**用法示例：**
```python
from datacommons_client import DataCommonsClient

client = DataCommonsClient()

# Get latest population for multiple entities
response = client.observation.fetch(
    variable_dcids=["Count_Person"],
    entity_dcids=["geoId/06", "geoId/48"],  # California and Texas
    date="latest"
)

# Get complete time series
response = client.observation.fetch(
    variable_dcids=["Count_Person"],
    entity_dcids=["country/USA"],
    date="all"
)

# Use relation expressions to query hierarchies
response = client.observation.fetch(
    variable_dcids=["Count_Person"],
    entity_expression="geoId/06<-containedInPlace+{typeOf:County}",
    date="2020"
)
```

### 2. fetch_available_statistical_variables()

发现哪些统计变量包含给定实体的数据。

**输入：** 仅实体 DCID
**输出：** 按实体组织的可用变量的字典

**用法示例：**
<<<代码块_1>>>

### 3. fetch_observations_by_entity_dcid()

通过 DCID 定位特定实体的显式方法（功能上相当于使用实体_dcids 的 `fetch()`）。

### 4. fetch_observations_by_entity_type()

检索按父级和类型分组的多个实体的观测值 - 对于查询一个区域中的所有国家或一个州内的所有县很有用。

**参数：**
- `parent_entity`：父实体 DCID
- `entity_type`：子实体的类型
- `variable_dcids`：要查询的统计变量
- `date`：时间规范
- `select` 和过滤器选项

**用法示例：**
<<<代码块_2>>>

## 响应对象方法

所有响应对象都支持：
- `to_json()`：格式为 JSON 字符串
- `to_dict()`：以字典形式返回
- `get_data_by_entity()`：按实体而不是变量重新组织
- `to_observations_as_records()`：展平为单独的记录

## 常见用例

### 用例 1：查询前检查数据可用性

使用 `select=["entity", "variable"]` 确认实体具有观察结果，而不检索实际数据：
<<<代码块_3>>>

### 用例 2：访问完整时间序列

请求 `date="all"` 获取完整的历史观察结果以进行趋势分析：
<<<代码块_4>>>

### 用例 3：按数据源过滤

指定 `filter_facet_domains` 从特定源检索数据以确保一致性：
<<<代码块_5>>>

### 用例 4：查询层次关系

使用关系表达式获取相关实体的观察结果：
<<<代码块_6>>>

## 与熊猫一起工作

该 API 与 Pandas 无缝集成。安装 Pandas 支持：
```bash
pip install "datacommons-client[Pandas]"
```

响应对象可以转换为 DataFrame 进行分析：
```python
response = client.observation.fetch(
    variable_dcids=["Count_Person"],
    entity_dcids=["geoId/06", "geoId/48"],
    date="all"
)

# Convert to DataFrame
df = response.to_observations_as_records()
# Returns DataFrame with columns: date, entity, variable, value
```

## 重要提示

- **facet** 代表数据源并包括来源元数据
- **orderedFacets** 按可靠性/新近度排序
- 使用关系表达式进行复杂的图形查询
- `fetch()` 方法是最灵活的 - 将其用于大多数查询