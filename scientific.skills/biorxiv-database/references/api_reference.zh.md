<!-- 此文件由机器翻译自 api_reference.md -->

# bioRxiv API 参考

## 概述

bioRxiv API 提供对来自 bioRxiv 服务器的预印本元数据的编程访问。该 API 返回 JSON 格式的数据，其中包含有关生命科学预印本的全面元数据。

## 基本网址

```
https://api.biorxiv.org
```

## 速率限制

尊重 API：
- 在请求之间添加延迟（建议至少 0.5 秒）
- 使用适当的用户代理标头
- 尽可能缓存结果

## API 端点

### 1. 按日期范围划分的详细信息

检索在特定日期范围内发布的预印本。

**端点：**
<<<代码块_1>>>

**参数：**
- `start_date`：开始日期，格式为 YYYY-MM-DD
- `end_date`：结束日期（YYYY-MM-DD 格式）
- `category`（可选）：按主题类别过滤

**示例：**
<<<代码块_2>>>

**回应：**
<<<代码块_3>>>

### 2. DOI 的详细信息

检索 DOI 特定预印本的详细信息。

**端点：**
<<<代码块_4>>>

**参数：**
- `doi`：预印本的 DOI（例如，`10.1101/2024.01.15.123456`）

**示例：**
<<<代码块_5>>>

### 3. 按时间间隔发布的出版物

从某个时间间隔检索最近的出版物。

**端点：**
<<<代码块_6>>>

**参数：**
- `interval`：返回搜索的天数（例如，过去 24 小时的 `1`）
- `cursor`：分页光标（第一页为 0，后续页递增 100）
- `format`：响应格式（`json` 或 `xml`）

**示例：**
```
GET https://api.biorxiv.org/pubs/biorxiv/1/0/json
```

**响应包括分页：**
```json
{
  "messages": [
    {
      "status": "ok",
      "count": 100,
      "total": 250,
      "cursor": 100
    }
  ],
  "collection": [...]
}
```

## 有效类别

bioRxiv 将预印本分为以下几类：

- `animal-behavior-and-cognition`
- `biochemistry`
- `bioengineering`
- `bioinformatics`
- `biophysics`
- `cancer-biology`
- `cell-biology`
- `clinical-trials`
- `developmental-biology`
- `ecology`
- `epidemiology`
- `evolutionary-biology`
- `genetics`
- `genomics`
- `immunology`
- `microbiology`
- `molecular-biology`
- `neuroscience`
- `paleontology`
- `pathology`
- `pharmacology-and-toxicology`
- `physiology`
- `plant-biology`
- `scientific-communication-and-education`
- `synthetic-biology`
- `systems-biology`
- `zoology`

## 论文元数据字段

`collection` 数组中的每篇论文包含：

|领域 |描述 |类型 |
|--------|-------------|------|
| `doi` |数字对象标识符|字符串|
| `title` |论文标题 |字符串|
| `authors` |以逗号分隔的作者列表 |字符串|
| `author_corresponding` |通讯作者姓名 |字符串|
| `author_corresponding_institution` |通讯作者单位 |字符串|
| `date` |出版日期 (YYYY-MM-DD) |字符串|
| `version` |版本号 |字符串|
| `type` |提交类型（例如“新结果”）|字符串|
| `license` |许可证类型（例如“cc_by”）|字符串|
| `category` |学科类别 |字符串|
| `jatsxml` | JATS XML 的 URL |字符串|
| `abstract` |论文摘要|字符串|
| `published` |期刊出版信息（如果已出版）|字符串|

## 下载全文

### PDF下载

PDF可以直接下载（不通过API）：

```
https://www.biorxiv.org/content/{doi}v{version}.full.pdf
```

示例：
```
https://www.biorxiv.org/content/10.1101/2024.01.15.123456v1.full.pdf
```

### HTML 版本

```
https://www.biorxiv.org/content/{doi}v{version}
```

### JATS XML

完整的结构化 XML 可通过 API 响应中的 `jatsxml` 字段获得。

## 常见搜索模式

### 作者搜索

1. 获取日期范围内的论文
2. 按作者姓名过滤（`authors` 字段中不区分大小写的子字符串匹配）

### 关键字搜索

1. 获取日期范围内的论文（可选择按类别过滤）
2. 在标题、摘要或这两个字段中搜索
3.过滤含有关键词的论文（不区分大小写）

### 最近的论文（按类别）

1.使用`/pubs/biorxiv/{interval}/0/json`端点
2. 如果需要按类别过滤

## 错误处理

常见的HTTP状态码：
- `200`：成功
- `404`：找不到资源
- `500`：服务器错误

始终检查响应中的 `messages` 数组：
```json
{
  "messages": [
    {
      "status": "ok",
      "count": 100
    }
  ]
}
```

## 最佳实践

1. **缓存结果**：存储检索到的论文，避免重复的API调用
2. **使用适当的日期范围**：日期范围越小返回速度越快
3. **按类别过滤**：减少数据传输和处理时间
4. **批处理**：下载多个PDF时，在请求之间添加延迟
5. **错误处理**：始终检查响应状态并优雅地处理错误
6. **版本跟踪**：注意论文可以有多个版本

## Python 用法示例

```python
from biorxiv_search import BioRxivSearcher

searcher = BioRxivSearcher(verbose=True)

# Search by keywords
papers = searcher.search_by_keywords(
    keywords=["CRISPR", "gene editing"],
    start_date="2024-01-01",
    end_date="2024-12-31",
    category="genomics"
)

# Search by author
papers = searcher.search_by_author(
    author_name="Smith",
    start_date="2023-01-01",
    end_date="2024-12-31"
)

# Get specific paper
paper = searcher.get_paper_details("10.1101/2024.01.15.123456")

# Download PDF
searcher.download_pdf("10.1101/2024.01.15.123456", "paper.pdf")
```

## 外部资源

- bioRxiv 主页：https://www.biorxiv.org/
- API 文档：https://api.biorxiv.org/
- JATS XML 规范：https://jats.nlm.nih.gov/