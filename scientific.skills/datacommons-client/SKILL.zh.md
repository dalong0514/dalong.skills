<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：datacommons-客户端
描述：与 Data Commons 合作，该平台提供对全球来源公共统计数据的编程访问。在处理人口数据、经济指标、健康统计数据、环境数据或通过 Data Commons 提供的任何公共数据集时，请使用此技能。适用于查询人口统计、GDP数据、失业率、疾病流行情况、地理实体解析以及探索统计实体之间的关系。
---

# 数据共享客户端

## 概述

提供对 Data Commons Python API v2 的全面访问，用于查询统计观察结果、探索知识图和解析实体标识符。数据共享将来自人口普查局、卫生组织、环境机构和其他权威来源的数据汇总到统一的知识图谱中。

## 安装

安装支持 Pandas 的 Data Commons Python 客户端：

```bash
uv pip install "datacommons-client[Pandas]"
```

对于没有 Pandas 的基本用法：
<<<代码块_1>>>

## 核心能力

Data Commons API 包含三个主要端点，每个端点均在专用参考文件中详细介绍：

### 1.观察端点-统计数据查询

查询实体的时间序列统计数据。请参阅 `references/observation.md` 以获取全面的文档。

**主要用例：**
- 检索人口、经济、健康或环境统计数据
- 访问历史时间序列数据进行趋势分析
- 查询层次结构数据（一个州的所有县、一个地区的所有国家）
- 比较多个实体的统计数据
- 按数据源过滤以保持一致性

**常见模式：**
<<<代码块_2>>>

### 2.节点端点-知识图谱探索

探索知识图谱中的实体关系和属性。请参阅 `references/node.md` 以获取全面的文档。

**主要用例：**
- 发现实体的可用属性
- 导航地理层次结构（父/子关系）
- 检索实体名称和元数据
- 探索实体之间的联系
- 列出图中的所有实体类型

**常见模式：**
<<<代码块_3>>>

### 3.解析端点-实体识别

将实体名称、坐标或外部 ID 转换为数据共享 ID (DCID)。请参阅 `references/resolve.md` 以获取全面的文档。

**主要用例：**
- 将地名转换为DCID以供查询
- 将坐标解析为地点
- 将维基数据 ID 映射到数据共享实体
- 处理不明确的实体名称

**常见模式：**
<<<代码块_4>>>

## 典型工作流程

大多数数据共享查询都遵循以下模式：

1. **解析实体**（如果以名称开头）：
   <<<代码块_5>>>

2. **发现可用变量**（可选）：
   <<<代码块_6>>>

3. **查询统计数据**：
   ```python
   response = client.observation.fetch(
       variable_dcids=["Count_Person", "UnemploymentRate_Person"],
       entity_dcids=dcids,
       date="latest"
   )
   ```

4. **处理结果**：
   ```python
   # As dictionary
   data = response.to_dict()

   # As Pandas DataFrame
   df = response.to_observations_as_records()
   ```

## 查找统计变量

统计变量在 Data Commons 中使用特定的命名模式：

**常见变量模式：**
- `Count_Person` - 总人口
- `Count_Person_Female` - 女性人口
- `UnemploymentRate_Person` - 失业率
- `Median_Income_Household` - 家庭收入中位数
- `Count_Death` - 死亡人数
- `Median_Age_Person` - 中位年龄

**发现方法：**
```python
# Check what variables are available for an entity
available = client.observation.fetch_available_statistical_variables(
    entity_dcids=["geoId/06"]
)

# Or explore via the web interface
# https://datacommons.org/tools/statvar
```

## 与熊猫一起工作

所有观察响应都与 Pandas 集成：

```python
response = client.observation.fetch(
    variable_dcids=["Count_Person"],
    entity_dcids=["geoId/06", "geoId/48"],
    date="all"
)

# Convert to DataFrame
df = response.to_observations_as_records()
# Columns: date, entity, variable, value

# Reshape for analysis
pivot = df.pivot_table(
    values='value',
    index='date',
    columns='entity'
)
```

## API 认证

**对于 datacommons.org（默认）：**
- 需要 API 密钥
- 通过环境变量设置：`export DC_API_KEY="your_key"`
- 或者初始化时传递：`client = DataCommonsClient(api_key="your_key")`
- 请求密钥：https://apikeys.datacommons.org/

**对于自定义 Data Commons 实例：**
- 无需 API 密钥
- 指定自定义端点：`client = DataCommonsClient(url="https://custom.datacommons.org")`

## 参考文档

每个端点的综合文档可在 `references/` 目录中找到：

- **`references/observation.md`**：完整的观察 API 文档，包含所有方法、参数、响应格式和常见用例
- **`references/node.md`**：用于图形探索、属性查询和层次结构导航的完整 Node API 文档
- **`references/resolve.md`**：用于实体识别和 DCID 解析的完整 Resolve API 文档
- **`references/getting_started.md`**：包含端到端示例和常见模式的快速入门指南

## 其他资源

- **官方文档**：https://docs.datacommons.org/api/python/v2/
- **统计变量浏览器**：https://datacommons.org/tools/statvar
- **数据共享浏览器**：https://datacommons.org/browser/
- **GitHub 存储库**：https://github.com/datacommonsorg/api-python
## 有效使用技巧

1. **始终从解析开始**：在查询数据之前将名称转换为DCID
2. **对层次结构使用关系表达式**：一次查询所有子项而不是单独查询
3. **首先检查数据可用性**：使用 `fetch_available_statistical_variables()` 查看可查询的内容
4. **利用 Pandas 集成**：将响应转换为 DataFrame 进行分析
5. **缓存解析**：如果重复查询相同的实体，则存储名称→DCID映射
6. **按facet过滤以确保一致性**：使用`filter_facet_domains`确保数据来自同一来源
7. **阅读参考文档**：每个端点在 `references/` 目录中都有大量文档