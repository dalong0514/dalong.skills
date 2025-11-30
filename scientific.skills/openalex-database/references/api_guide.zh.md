<!-- 此文件由机器翻译自 api_guide.md -->

# OpenAlex API 完整指南

## 基础信息

**基本网址：** `https://api.openalex.org`
**身份验证：** 不需要
**速率限制：**
- 默认：1 个请求/秒，100k 个请求/天
- 礼貌池（带电子邮件）：10 个请求/秒，100k 个请求/天

## 关键最佳实践

### ✅ DO：使用 `?sample` 参数进行随机采样
```
https://api.openalex.org/works?sample=20&seed=123
```
对于大型样本 (10k+)，请使用多个种子并进行重复数据删除。

### ❌不要：使用随机页码进行采样
不正确：`?page=5`、`?page=17` - 这会使结果产生偏差！

### ✅ DO：使用两步查找进行实体过滤
<<<代码块_1>>>

### ❌不要：直接按实体名称过滤
错误：`/works?filter=author_name:Einstein` - 名称不明确！

### ✅ DO：使用最大页面大小进行批量提取
<<<代码块_2>>>
这比默认值 (25) 快 8 倍。

### ❌不要：使用默认页面尺寸
默认每页仅 25 个结果。

### ✅ DO：使用 OR 过滤器（管道 |）进行批量查找
<<<代码块_3>>>
每个过滤器最多 50 个值。

### ❌不要：对列表进行连续的 API 调用
当您可以批量调用时，进行 100 个单独的调用效率很低。

### ✅ 应：实施指数退避重试
<<<代码块_4>>>

### ✅ 要做：添加电子邮件以获得 10 倍速率限制提升
<<<代码块_5>>>
从 1 请求/秒增加到 10 请求/秒。

## 实体端点

- `/works` - 2.4 亿+ 学术文档
- `/authors` - 研究人员简介
- `/sources` - 期刊、存储库、会议
- `/institutions` - 大学、研究机构
- `/topics` - 主题分类（3 级层次结构）
- `/publishers` - 出版组织
- `/funders` - 资助机构
- `/text` - 使用主题/关键字标记您自己的文本 (POST)

## 基本查询参数

|参数|描述 |示例|
|------------|-------------|---------|
| `filter=` |筛选结果 | `?filter=publication_year:2020` |
| `search=` |全文检索 | `?search=machine+learning` |
| `sort=` |对结果排序 | `?sort=cited_by_count:desc` |
| `per-page=` |每页结果（最多 200 个）| `?per-page=200` |
| `page=` |页码 | `?page=2` |
| `sample=` |随机结果 | `?sample=50&seed=42` |
| `select=` |限制字段​​ | `?select=id,title` |
| `group_by=` |按领域聚合 | `?group_by=publication_year` |
| `mailto=` |礼貌泳池的电子邮件 | `?mailto=you@example.edu` |

## 过滤器语法

### 基本过滤
<<<代码块_6>>>

### 比较运算符
```
Greater than:      ?filter=cited_by_count:>100
Less than:         ?filter=publication_year:<2020
Range:             ?filter=publication_year:2020-2023
```

### 同一属性中有多个值
```
Repeat filter:     ?filter=institutions.country_code:us,institutions.country_code:gb
Use + symbol:      ?filter=institutions.country_code:us+gb
```
两者的意思都是：“与美国作者和英国作者合作”

### OR 查询
```
Any of these:      ?filter=institutions.country_code:us|gb|ca
Batch IDs:         ?filter=doi:10.1/abc|10.2/def
```
最多 50 个管道值。

## 常见查询模式

### 获取随机样本
```bash
# Small sample
https://api.openalex.org/works?sample=20&seed=42

# Large sample (10k+) - make multiple requests
https://api.openalex.org/works?sample=1000&seed=1
https://api.openalex.org/works?sample=1000&seed=2
# Then deduplicate by ID
```

### 搜索作品
```bash
# Simple search
https://api.openalex.org/works?search=machine+learning

# Search specific field
https://api.openalex.org/works?filter=title.search:CRISPR

# Search + filter
https://api.openalex.org/works?search=climate&filter=publication_year:2023
```

### 查找作者作品（两步）
```bash
# Step 1: Get author ID
https://api.openalex.org/authors?search=Heather+Piwowar
# Returns: "id": "https://openalex.org/A5023888391"

# Step 2: Get their works
https://api.openalex.org/works?filter=authorships.author.id:A5023888391
```

### 按机构查找作品（两步）
```bash
# Step 1: Get institution ID
https://api.openalex.org/institutions?search=MIT
# Returns: "id": "https://openalex.org/I136199984"

# Step 2: Get their works
https://api.openalex.org/works?filter=authorships.institutions.id:I136199984
```

### 最近被引用次数最多的论文
```bash
https://api.openalex.org/works?filter=publication_year:>2020&sort=cited_by_count:desc&per-page=200
```

### 开放获取作品
```bash
# All OA
https://api.openalex.org/works?filter=is_oa:true

# Gold OA only
https://api.openalex.org/works?filter=open_access.oa_status:gold
```

### 多重标准
```bash
# Recent OA works about COVID from top institutions
https://api.openalex.org/works?filter=publication_year:2022,is_oa:true,title.search:covid,authorships.institutions.id:I136199984|I27837315
```

### 批量 DOI 查询
```bash
# Get specific works by DOI (up to 50 per request)
https://api.openalex.org/works?filter=doi:https://doi.org/10.1371/journal.pone.0266781|https://doi.org/10.1371/journal.pone.0267149&per-page=50
```

### 汇总数据
```bash
# Top topics
https://api.openalex.org/works?group_by=topics.id

# Papers per year
https://api.openalex.org/works?group_by=publication_year

# Most prolific institutions
https://api.openalex.org/works?group_by=authorships.institutions.id
```

### 分页
```bash
# First page
https://api.openalex.org/works?filter=publication_year:2023&per-page=200

# Next pages
https://api.openalex.org/works?filter=publication_year:2023&per-page=200&page=2
```

## 响应结构

### 列出端点
```json
{
  "meta": {
    "count": 240523418,
    "db_response_time_ms": 42,
    "page": 1,
    "per_page": 25
  },
  "results": [
    { /* entity object */ }
  ]
}
```

### 单一实体
```
https://api.openalex.org/works/W2741809807
→ Returns Work object directly (no meta/results wrapper)
```

### 分组依据
```json
{
  "meta": { "count": 100 },
  "group_by": [
    {
      "key": "https://openalex.org/T10001",
      "key_display_name": "Artificial Intelligence",
      "count": 15234
    }
  ]
}
```

## 工作过滤器（最常见）

|过滤|描述 |示例|
|--------|-------------|---------|
| `authorships.author.id` |作者的 OpenAlex ID | `A5023888391` |
| `authorships.institutions.id` |机构 ID | `I136199984` |
| `cited_by_count` |引用计数 | `>100` |
| `is_oa` |是开放获取 | `true/false` |
| `publication_year` |出版年份 | `2020`、`>2020`、`2018-2022` |
| `primary_location.source.id` |来源（期刊）ID | `S137773608` |
| `topics.id` |主题 ID | `T10001` |
| `type` |文件类型 | `article`、`book`、`dataset` |
| `has_doi` |拥有 DOI | `true/false` |
| `has_fulltext` |有全文 | `true/false` |

## 作者过滤器

|过滤|描述 |
|--------|-------------|
| `last_known_institution.id` |当前/最后一个机构 |
| `works_count` |作品数量 |
| `cited_by_count` |总引用次数 |
| `orcid` | ORCID 标识符 |

## 外部 ID 支持

### 作品
```
DOI:  /works/https://doi.org/10.7717/peerj.4375
PMID: /works/pmid:29844763
```

### 作者
```
ORCID: /authors/https://orcid.org/0000-0003-1613-5981
```

### 机构
```
ROR: /institutions/https://ror.org/02y3ad647
```

### 来源
```
ISSN: /sources/issn:0028-0836
```

## 性能提示

1. **使用最大页面大小**：`?per-page=200`（调用次数减少 8 倍）
2. **批量ID查找**：使用管道运算符最多查找50个ID
3. **仅选择需要的字段**：`?select=id,title,publication_year`
4. **使用并发请求**：具有速率限制（对于电子邮件，每秒 10 个请求）
5. **添加电子邮件**：`?mailto=you@example.edu` 速度提升 10 倍

## 错误处理

### HTTP 状态代码
- `200` - 成功
- `400` - 错误请求（检查过滤器语法）
- `403` - 超出速率限制（实施退避）
- `404` - 实体不存在
- `500` - 服务器错误（退避重试）

### 指数退避
```python
def fetch_with_retry(url, max_retries=5):
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                return response.json()
            elif response.status_code in [403, 500, 502, 503, 504]:
                wait_time = 2 ** attempt
                time.sleep(wait_time)
            else:
                response.raise_for_status()
        except requests.exceptions.Timeout:
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
            else:
                raise
    raise Exception(f"Failed after {max_retries} retries")
```

## 速率限制

### 没有电子邮件（默认池）
- 1 个请求/秒
- 100,000 个请求/天

### 通过电子邮件（礼貌池）
- 10 个请求/秒
- 100,000 个请求/天
- **始终用于生产**

### 并发请求策略
1. 跟踪全球每秒请求数
2.跨线程使用信号量或速率限制器
3.监控403响应
4.如果达到限制就退缩

## 要避免的常见错误

1. ❌ 使用页码进行采样 → ✅ 使用 `?sample=`
2. ❌ 按实体名称过滤 → ✅ 先获取ID
3. ❌ 默认页面大小 → ✅ 使用 `per-page=200`
4. ❌ 顺序 ID 查找 → ✅ 使用管道运算符进行批处理
5. ❌ 无错误处理 → ✅ 实施带退避的重试
6. ❌ 忽略速率限制 → ✅ 全局速率限制
7. ❌ 不包括电子邮件 → ✅ 添加 `mailto=`
8. ❌ 获取所有字段 → ✅ 使用 `select=`

## 其他资源

- 完整文档：https://docs.openalex.org
- API 概述：https://docs.openalex.org/how-to-use-the-api/api-overview
- 实体模式：https://docs.openalex.org/api-entities
- 帮助：https://openalex.org/help
- 用户组：https://groups.google.com/g/openalex-users