<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：openalex-数据库
描述：使用 OpenAlex 数据库查询和分析学术文献。在搜索学术论文、分析研究趋势、查找作者或机构的作品、跟踪引用、发现开放获取出版物或对 2.4 亿多篇学术作品进行文献计量分析时，应该使用此技能。用于文献检索、研究成果分析、引文分析和学术数据库查询。
---

# OpenAlex 数据库

## 概述

OpenAlex 是一个全面的开放目录，包含 2.4 亿多篇学术著作、作者、机构、主题、来源、出版商和资助者。该技能提供了用于查询 OpenAlex API 来搜索文献、分析研究成果、跟踪引用和进行文献计量研究的工具和工作流程。

## 快速入门

### 基本设置

始终使用电子邮件地址初始化客户端以访问礼貌池（10 倍速率限制提升）：

```python
from scripts.openalex_client import OpenAlexClient

client = OpenAlexClient(email="your-email@example.edu")
```

### 安装要求

使用 uv 安装所需的包：

<<<代码块_1>>>

不需要 API 密钥 - OpenAlex 是完全开放的。

## 核心能力

### 1. 搜索论文

**用途**：按标题、摘要或主题查找论文

<<<代码块_2>>>

### 2. 按作者查找作品

**用途**：获取特定研究人员的所有出版物

使用两步模式（实体名称 → ID → 作品）：

<<<代码块_3>>>

**手动两步方法**：
<<<代码块_4>>>

### 3.从机构查找作品

**用途**：分析大学或组织的研究成果

<<<代码块_5>>>

### 4.高被引论文

**用途**：寻找某个领域有影响力的论文

<<<代码块_6>>>

### 5. 开放获取论文

**用途**：查找免费可用的研究

```python
from scripts.query_helpers import get_open_access_papers

papers = get_open_access_papers(
    search_term="climate change",
    client=client,
    oa_status="any",  # or "gold", "green", "hybrid", "bronze"
    limit=200
)
```

### 6. 出版趋势分析

**用途**：跟踪一段时间内的研究成果

```python
from scripts.query_helpers import get_publication_trends

trends = get_publication_trends(
    search_term="artificial intelligence",
    filter_params={"is_oa": "true"},
    client=client
)

# Sort and display
for trend in sorted(trends, key=lambda x: x['key'])[-10:]:
    print(f"{trend['key']}: {trend['count']} publications")
```

### 7. 研究成果分析

**用途**：作者或机构研究的综合分析

```python
from scripts.query_helpers import analyze_research_output

analysis = analyze_research_output(
    entity_type='institution',  # or 'author'
    entity_name='MIT',
    client=client,
    years='>2020'
)

print(f"Total works: {analysis['total_works']}")
print(f"Open access: {analysis['open_access_percentage']}%")
print(f"Top topics: {analysis['top_topics'][:5]}")
```

### 8. 批量查找

**用途**：高效获取多个 DOI、ORCID 或 ID 的信息

```python
dois = [
    "https://doi.org/10.1038/s41586-021-03819-2",
    "https://doi.org/10.1126/science.abc1234",
    # ... up to 50 DOIs
]

works = client.batch_lookup(
    entity_type='works',
    ids=dois,
    id_field='doi'
)
```

### 9. 随机抽样

**用途**：获取代表性样品进行分析

```python
# Small sample
works = client.sample_works(
    sample_size=100,
    seed=42,  # For reproducibility
    filter_params={"publication_year": "2023"}
)

# Large sample (>10k) - automatically handles multiple requests
works = client.sample_works(
    sample_size=25000,
    seed=42,
    filter_params={"is_oa": "true"}
)
```

### 10. 引文分析

**用途**：查找引用特定作品的论文

```python
# Get the work
work = client.get_entity('works', 'https://doi.org/10.1038/s41586-021-03819-2')

# Get citing papers using cited_by_api_url
import requests
citing_response = requests.get(
    work['cited_by_api_url'],
    params={'mailto': client.email, 'per-page': 200}
)
citing_works = citing_response.json()['results']
```

### 11.主题和主题分析

**用途**：了解研究重点领域

```python
# Get top topics for an institution
topics = client.group_by(
    entity_type='works',
    group_field='topics.id',
    filter_params={
        "authorships.institutions.id": "I136199984",  # MIT
        "publication_year": ">2020"
    }
)

for topic in topics[:10]:
    print(f"{topic['key_display_name']}: {topic['count']} works")
```

### 12.大规模数据提取

**用途**：下载大型数据集进行分析

```python
# Paginate through all results
all_papers = client.paginate_all(
    endpoint='/works',
    params={
        'search': 'synthetic biology',
        'filter': 'publication_year:2020-2024'
    },
    max_results=10000
)

# Export to CSV
import csv
with open('papers.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['Title', 'Year', 'Citations', 'DOI', 'OA Status'])

    for paper in all_papers:
        writer.writerow([
            paper.get('title', 'N/A'),
            paper.get('publication_year', 'N/A'),
            paper.get('cited_by_count', 0),
            paper.get('doi', 'N/A'),
            paper.get('open_access', {}).get('oa_status', 'closed')
        ])
```

## 关键最佳实践

### 始终使用电子邮件以保持礼貌
添加电子邮件以获得 10 倍速率限制（1 请求/秒 → 10 请求/秒）：
```python
client = OpenAlexClient(email="your-email@example.edu")
```

### 使用两步模式进行实体查找
切勿直接按实体名称过滤 - 始终首先获取 ID：
```python
# ✅ Correct
# 1. Search for entity → get ID
# 2. Filter by ID

# ❌ Wrong
# filter=author_name:Einstein  # This doesn't work!
```

### 使用最大页面大小
始终使用 `per-page=200` 进行高效的数据检索：
```python
results = client.search_works(search="topic", per_page=200)
```

### 批量多个ID
对多个ID而不是单个请求使用batch_lookup()：
```python
# ✅ Correct - 1 request for 50 DOIs
works = client.batch_lookup('works', doi_list, 'doi')

# ❌ Wrong - 50 separate requests
for doi in doi_list:
    work = client.get_entity('works', doi)
```

### 对随机数据使用样本参数
将 `sample_works()` 与种子一起使用以实现可重复的随机采样：
```python
# ✅ Correct
works = client.sample_works(sample_size=100, seed=42)

# ❌ Wrong - random page numbers bias results
# Using random page numbers doesn't give true random sample
```

### 仅选择需要的字段
通过选择特定字段来减少响应大小：
```python
results = client.search_works(
    search="topic",
    select=['id', 'title', 'publication_year', 'cited_by_count']
)
```

## 常见过滤器模式

### 日期范围
```python
# Single year
filter_params={"publication_year": "2023"}

# After year
filter_params={"publication_year": ">2020"}

# Range
filter_params={"publication_year": "2020-2024"}
```

### 多个过滤器（AND）
```python
# All conditions must match
filter_params={
    "publication_year": ">2020",
    "is_oa": "true",
    "cited_by_count": ">100"
}
```

### 多个值（或）
```python
# Any institution matches
filter_params={
    "authorships.institutions.id": "I136199984|I27837315"  # MIT or Harvard
}
```

### 协作（以及属性内）
```python
# Papers with authors from BOTH institutions
filter_params={
    "authorships.institutions.id": "I136199984+I27837315"  # MIT AND Harvard
}
```

### 否定
```python
# Exclude type
filter_params={
    "type": "!paratext"
}
```

## 实体类型

OpenAlex 提供以下实体类型：
- **作品** - 学术文档（文章、书籍、数据集）
- **作者** - 身份明确的研究人员
- **机构** - 大学和研究机构
- **来源** - 期刊、知识库、会议
- **主题** - 主题分类
- **出版商** - 出版组织
- **资助者** - 资助机构

使用一致的模式访问任何实体类型：
```python
client.search_works(...)
client.get_entity('authors', author_id)
client.group_by('works', 'topics.id', filter_params={...})
```

## 外部 ID

直接使用外部标识符：
```python
# DOI for works
work = client.get_entity('works', 'https://doi.org/10.7717/peerj.4375')

# ORCID for authors
author = client.get_entity('authors', 'https://orcid.org/0000-0003-1613-5981')

# ROR for institutions
institution = client.get_entity('institutions', 'https://ror.org/02y3ad647')

# ISSN for sources
source = client.get_entity('sources', 'issn:0028-0836')
```

## 参考文档

### 详细 API 参考
请参阅`references/api_guide.md`了解：
- 完整的过滤器语法
- 所有可用端点
- 响应结构
- 错误处理
- 性能优化
- 速率限制细节

### 常见查询示例
请参阅 `references/common_queries.md` 了解：
- 完整的工作示例
- 真实世界的用例
- 复杂的查询模式
- 数据导出工作流程
- 多步骤分析程序

## 脚本
### openalex_client.py
主要 API 客户端：
- 自动速率限制
- 指数退避重试逻辑
- 分页支持
- 批量操作
- 错误处理

用于完全控制的直接 API 访问。

### query_helpers.py
用于常见操作的高级辅助函数：
- `find_author_works()` - 按作者获取论文
- `find_institution_works()` - 从机构获取论文
- `find_highly_cited_recent_papers()` - 获取有影响力的论文
- `get_open_access_papers()` - 查找 OA 出版物
- `get_publication_trends()` - 分析一段时间内的趋势
- `analyze_research_output()` - 综合分析

用于通过简化的界面进行常见的研究查询。

## 故障排除

### 速率限制
如果遇到 403 错误：
1. 确保电子邮件已添加到请求中
2. 验证不超过 10 请求/秒
3.客户端自动实现指数退避

### 空结果
如果搜索没有返回结果：
1. 检查过滤器语法（参见`references/api_guide.md`）
2. 使用两步模式进行实体查找（不按名称过滤）
3.验证实体ID的格式是否正确

### 超时错误
对于大型查询：
1. 使用 `per-page=200` 进行分页
2.使用`select=`限制返回字段
3. 如果需要，分成更小的查询

## 速率限制

- **默认**：1 个请求/秒，100k 个请求/天
- **礼貌池（带电子邮件）**：10 个请求/秒，100k 请求/天

通过向客户提供电子邮件，始终使用礼貌池进行生产工作流程。

## 注释

- 无需身份验证
- 所有数据都是开放且免费的
- 速率限制适用于全球，而不是每个 IP
- 如果需要基于 LLM 的分析，请使用 LitLLM 和 OpenRouter（不要直接使用 Perplexity API）
- 客户端自动处理分页、重试和速率限制