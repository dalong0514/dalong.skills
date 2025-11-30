<!-- 此文件由机器翻译自 api_reference.md -->

# PubMed 电子实用程序 API 参考

## 概述

NCBI 电子实用程序通过 REST API 提供对 PubMed 和其他 Entrez 数据库的编程访问。所有电子实用程序的基本 URL 是：

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
```

## API 密钥要求

自 2018 年 12 月 1 日起，NCBI 强制要求电子实用程序调用使用 API 密钥。 API 密钥将速率限制从 3 个请求/秒增加到 10 个请求/秒。要获取 API 密钥，请注册 NCBI 帐户并从您的帐户设置生成密钥。

使用 `&api_key` 参数在请求中包含 API 密钥：
<<<代码块_1>>>

## 速率限制

- **没有 API 密钥**：每秒 3 个请求
- **使用 API 密钥**：每秒 10 个请求
- 始终在请求中包含 User-Agent 标头

## 核心电子实用工具

### 1. ESearch - 查询数据库

**端点**：`esearch.fcgi`

**用途**：搜索 Entrez 数据库并检索 UID 列表（例如 PubMed 的 PMID）

**所需参数**：
- `db` - 要搜索的数据库（例如，pubmed、基因、蛋白质）
- `term` - 搜索查询

**可选参数**：
- `retmax` - 返回的最大记录数（默认值：20，最大值：10000）
- `retstart` - 要返回的第一条记录的索引（默认值：0）
- `usehistory=y` - 将大型结果集的结果存储在历史服务器上
- `retmode` - 返回格式（xml、json）
- `sort` - 排序顺序（相关性、pub_date、first_author、last_author、期刊）
- `field` - 将搜索限制为特定字段
- `datetype` - 用于过滤的日期类型（pdat 表示发布日期）
- `mindate` - 最短日期（YYYY/MM/DD 格式）
- `maxdate` - 最大日期（YYYY/MM/DD 格式）

**请求示例**：
<<<代码块_2>>>

**响应元素**：
- `Count` - 匹配查询的记录总数
- `RetMax` - 此响应中返回的记录数
- `RetStart` - 第一个返回记录的索引
- `IdList` - UID 列表 (PMID)
- `WebEnv` - 历史服务器环境字符串（当 usehistory=y 时）
- `QueryKey` - 历史服务器的查询键（当usehistory=y时）

### 2. EFetch - 下载记录

**端点**：`efetch.fcgi`

**目的**：从数据库中检索各种格式的完整记录

**所需参数**：
- `db` - 数据库名称
- `id` - 以逗号分隔的 UID 列表，或使用 ESearch 中的 WebEnv/query_key

**可选参数**：
- `rettype` - 记录类型（abstract、medline、xml、uilist）
- `retmode` - 返回模式（文本、xml）
- `retstart` - 起始记录索引
- `retmax` - 每个请求的最大记录数

**请求示例**：
<<<代码块_3>>>

**PubMed 的常见 rettype 值**：
- `abstract` - 摘要文本
- `medline` - 完整 MEDLINE 格式
- `xml` - PubMed XML 格式
- `uilist` - 仅 UID 列表

### 3. ESummary - 检索文档摘要

**端点**：`esummary.fcgi`

**用途**：获取 UID 列表的文档摘要 (DocSum)

**所需参数**：
- `db` - 数据库名称
- `id` - 逗号分隔的 UID 或 WebEnv/query_key

**可选参数**：
- `retmode` - 返回格式（xml、json）
- `version` - DocSum 版本（1.0 或 2.0，默认为 1.0）

**请求示例**：
<<<代码块_4>>>

**DocSum 字段**（因数据库而异，常见的 PubMed 字段）：
- 标题
- 作者
- 来源（期刊）
- 发布日期
- 卷、期、页数
- DOI
- PmcRefCount（PMC 中的引用）

### 4. EPost - 上传 UID

**端点**：`epost.fcgi`

**用途**：将UID列表上传到历史服务器以供后续请求使用

**所需参数**：
- `db` - 数据库名称
- `id` - 以逗号分隔的 UID 列表

**请求示例**：
<<<代码块_5>>>

**回应**：
返回WebEnv和QueryKey以供后续请求使用

### 5. ELink - 查找相关数据

**端点**：`elink.fcgi`

**目的**：查找同一数据库内或不同数据库内的相关记录

**所需参数**：
- `dbfrom` - 源数据库
- `db` - 目标数据库（可以与 dbfrom 相同）
- `id` - 来自源数据库的 UID

**可选参数**：
- `cmd` - 链接命令（neighbor、neighbor_history、prlinks、llinks 等）
- `linkname` - 要检索的特定链接类型
- `term` - 使用搜索查询过滤结果
- `holding` - 按图书馆馆藏过滤

**请求示例**：
<<<代码块_6>>>

**常用链接命令**：
- `neighbor` - 返回相关记录
- `neighbor_history` - 将相关记录发布到历史服务器
- `prlinks` - 返回提供商 URL
- `llinks` - 返回 LinkOut URL

### 6. EInfo - 数据库信息

**端点**：`einfo.fcgi`

**用途**：获取有关可用 Entrez 数据库或特定数据库字段的信息

**参数**：
- `db` - 数据库名称（可选；省略列出所有数据库）
- `retmode` - 返回格式（xml、json）

**请求示例**：
```
einfo.fcgi?db=pubmed&retmode=json&api_key=YOUR_API_KEY
```

**回报**：
- 数据库描述
- 记录数
- 最后更新日期
- 带有描述的可用搜索字段

### 7. EGQuery - 全局查询

**端点**：`egquery.fcgi`

**目的**：所有 Entrez 数据库中的搜索词计数

**所需参数**：
- `term` - 搜索查询

**请求示例**：
```
egquery.fcgi?term=cancer&api_key=YOUR_API_KEY
```

### 8.ESpell - 拼写建议

**端点**：`espell.fcgi`

**目的**：获取查询的拼写建议

**所需参数**：
- `db` - 数据库名称
- `term` - 可能存在拼写错误的搜索词

**请求示例**：
```
espell.fcgi?db=pubmed&term=cancre&api_key=YOUR_API_KEY
```

### 9. ECitMatch - 引文匹配

**端点**：`ecitmatch.cgi`

**目的**：使用期刊、年份、卷、页数、作者信息搜索 PubMed 引文

**请求格式**：带有引用字符串的 POST 请求

**引文字符串格式**：
```
journal|year|volume|page|author|key|
```

**示例**：
```
Science|2008|320|5880|1185|key1|
Nature|2010|463|7279|318|key2|
```

**速率限制**：每秒 3 个请求，需要 User-Agent 标头

## 最佳实践

### 使用历史服务器处理大型结果集

对于返回超过 500 条记录的查询，请使用历史记录服务器：

1. **初始搜索历史记录**：
```
esearch.fcgi?db=pubmed&term=cancer&usehistory=y&retmode=json&api_key=YOUR_API_KEY
```

2. **批量检索记录**：
```
efetch.fcgi?db=pubmed&query_key=1&WebEnv=MCID_12345&retstart=0&retmax=500&rettype=xml&api_key=YOUR_API_KEY
efetch.fcgi?db=pubmed&query_key=1&WebEnv=MCID_12345&retstart=500&retmax=500&rettype=xml&api_key=YOUR_API_KEY
```

### 批量操作

在获取之前使用 EPost 上传大量 UID 列表：

```
# Step 1: Post UIDs
epost.fcgi?db=pubmed&id=123,456,789,...&api_key=YOUR_API_KEY

# Step 2: Fetch using WebEnv/query_key
efetch.fcgi?db=pubmed&query_key=1&WebEnv=MCID_12345&rettype=xml&api_key=YOUR_API_KEY
```

### 错误处理

常见的HTTP状态码：
- `200` - 成功
- `400` - 错误请求（检查参数）
- `414` - URI 太长（使用 POST 或历史服务器）
- `429` - 超出速率限制

### 缓存

实施本地缓存以：
- 减少冗余的API调用
- 保持在费率限制范围内
- 改善响应时间
- 尊重 NCBI 资源

## 响应格式

### XML（默认）

最详细的格式与完整的结构化数据。每个数据库都有自己的DTD（文档类型定义）。

### JSON

适用于大多数具有 `retmode=json` 的实用程序。在现代应用程序中更容易解析。

### 文字

纯文本格式，适用于摘要和简单的数据检索。

## 支持和资源

- **API 文档**：https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **邮件列表**：utilities-announce@ncbi.nlm.nih.gov
- **支持**：vog.hin.mln.ibcn@seitilitue
- **NLM 服务台**：1-888-FIND-NLM (1-888-346-3656)