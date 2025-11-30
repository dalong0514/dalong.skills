<!-- 此文件由机器翻译自 patentsearch_api.md -->

# PatentSearch API 参考

## 概述

PatentSearch API 是 USPTO 基于 ElasticSearch 的现代专利检索系统，于 2025 年 5 月取代了旧版 PatentsView API。它提供对 2025 年 6 月 30 日之前的专利数据的访问，并定期更新。

**基本网址：** `https://search.patentsview.org/api/v1/`

## 身份验证

所有 API 请求都需要使用请求标头中的 API 密钥进行身份验证：

```
X-Api-Key: YOUR_API_KEY
```

注册 API 密钥：https://account.uspto.gov/api-manager/

## 速率限制

- **每个 API 密钥每分钟 45 个请求**
- 超出速率限制会导致 HTTP 429 错误

## 可用端点

### 核心专利和出版物端点

- **`/patent`** - 一般专利数据（已授权专利）
- **`/publication`** - 授权前发布数据
- **`/publication/rel_app_text`** - 出版物的相关申请数据

### 实体端点

- **`/inventor`>** - 包含位置和性别代码字段的发明人信息
- **`/assignee`>** - 带有位置标识符的受让人详细信息
- **`/location`>** - 地理数据，包括纬度/经度坐标
- **`/attorney`** - 法定代表人信息

### 分类端点

- **`/cpc_subclass`>** - 子类级别的合作专利分类
- **`/cpc_at_issue`** - 截至专利发布日期的 CPC 分类
- **`/uspc`** - 美国专利分类数据
- **`/wipo`** - 世界知识产权组织分类
- **`/ipc`** - 国际专利分类

### 文本数据端点（测试版）

- **`/brief_summary_text`** - 专利摘要（已授权和预授权）
- **`/claims`** - 专利权利要求文本
- **`/drawing_description_text`** - 绘图描述
- **`/detail_description_text`** - 详细描述文本

*注意：文本端点处于测试阶段，数据主要来自 2023 年以后。历史回填正在进行中。*

### 支持端点

- **`/other_reference`** - 专利参考资料
- **`/related_document`** - 专利之间的交叉引用

## 查询参数

所有端点都支持四个主要参数：

### 1. 查询字符串 (`q`)

使用 JSON 查询对象过滤数据。 **必填参数。**

**查询运算符：**

- **平等**：`{"field": "value"}` 或 `{"field": {"_eq": "value"}}`
- **不等于**：`{"field": {"_neq": "value"}}`
- **比较**：`_gt`、`_gte`、`_lt`、`_lte`
- **字符串匹配**：
  - `_begins` - 开头
  - `_contains` - 子字符串匹配
- **全文搜索**（推荐用于文本字段）：
  - `_text_all` - 所有术语必须匹配
  - `_text_any` - 任何术语匹配
  - `_text_phrase` - 精确短语匹配
- **逻辑运算符**：`_and`、`_or`、`_not`
- **数组匹配**：使用数组进行 OR 条件

**示例：**

<<<代码块_1>>>

### 2. 字段列表 (`f`)

指定在响应中返回哪些字段。可选 - 每个端点都有默认字段。

**格式：** 字段名称的 JSON 数组

<<<代码块_2>>>

### 3.排序(`s`)

按指定字段对结果进行排序。选修的。

**格式：** 带有字段名称和方向的 JSON 数组

<<<代码块_3>>>

### 4. 选项 (`o`)

控制分页和其他设置。选修的。

**可用选项：**

- `page` - 页码（默认值：1）
- `per_page` - 每页记录数（默认值：100，最大：1,000）
- `pad_patent_id` - 用前导零填充专利 ID（默认值： false）
- `exclude_withdrawn` - 排除撤回的专利（默认值：true）

**格式：** JSON 对象

<<<代码块_4>>>

## 响应格式

所有响应都遵循以下结构：

<<<代码块_5>>>

- `error` - 指示是否发生错误的布尔值
- `count` - 当前响应中的记录数
- `total_hits` - 匹配记录总数
- 端点特定的数据数组（例如，`patents`、`inventors`）

## 完整的请求示例

### 使用卷曲

<<<代码块_6>>>

### 使用Python

```python
import requests

url = "https://search.patentsview.org/api/v1/patent"
headers = {
    "X-Api-Key": "YOUR_API_KEY",
    "Content-Type": "application/json"
}
data = {
    "q": {
        "_and": [
            {"patent_date": {"_gte": "2024-01-01"}},
            {"patent_abstract": {"_text_all": ["artificial", "intelligence"]}}
        ]
    },
    "f": ["patent_number", "patent_title", "patent_date", "assignee_organization"],
    "s": [{"patent_date": "desc"}],
    "o": {"per_page": 100}
}

response = requests.post(url, headers=headers, json=data)
results = response.json()
```

## 常用字段名称

### 专利端点字段

- `patent_number` - 专利号
- `patent_title` - 专利标题
- `patent_date` - 授予日期
- `patent_abstract` - 摘要文本
- `patent_type` - 专利类型
- `inventor_name` - 发明人姓名（数组）
- `assignee_organization` - 受让人公司名称（数组）
- `cpc_subclass_id` - CPC 分类代码
- `uspc_class` - 美国分类代码
- `cited_patent_number` - 对其他专利的引用
- `citedby_patent_number` - 引用本专利的专利

请参阅完整字段字典：https://search.patentsview.org/docs/

## 最佳实践

1. **对文本字段使用 `_text*` 运算符** - 比 `_contains` 或 `_begins` 性能更高
2. **仅请求所需的字段** - 减少响应大小并提高性能
3. **实现分页** - 高效处理大型结果集
4. **尊重率限制** - 针对 429 错误实施退避/重试逻辑
5. **缓存结果** - 减少冗余的API调用
6. **使用日期范围** - 缩小搜索范围以提高性能

## 错误处理

常见的HTTP状态码：

- **200** - 成功
- **400** - 错误请求（无效的查询语法）
- **401** - 未经授权（API 密钥丢失或无效）
- **429** - 请求过多（超出速率限制）
- **500** - 服务器错误

## 最近更新（2025 年 2 月）

- 数据更新至 2024 年 12 月 31 日
- 用于格式化专利 ID 的新 `pad_patent_id` 选项
- 新的`exclude_withdrawn`选项显示撤回的专利
- 文本端点继续测试回填

## 资源

- **官方文档**：https://search.patentsview.org/docs/
- **API 密钥注册**：https://account.uspto.gov/api-manager/
- **旧版 API 通知**：旧版 PatentsView API 已于 2025 年 5 月 1 日停用