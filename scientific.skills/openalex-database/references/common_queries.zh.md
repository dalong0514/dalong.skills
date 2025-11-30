<!-- 此文件由机器翻译自 common_queries.md -->

# 常见 OpenAlex 查询示例

本文档提供了使用 OpenAlex 进行常见研究查询的实际示例。

## 按作者查找论文

**用户查询**：“查找阿尔伯特·爱因斯坦的论文”

**方法**：两步模式
1.搜索作者获取ID
2.按作者ID过滤作品

**Python 示例**：
```python
from scripts.openalex_client import OpenAlexClient
from scripts.query_helpers import find_author_works

client = OpenAlexClient(email="your-email@example.edu")
works = find_author_works("Albert Einstein", client, limit=100)

for work in works:
    print(f"{work['title']} ({work['publication_year']})")
```

## 从机构查找论文

**用户查询**：“麻省理工学院去年发表了哪些论文？”

**方法**：带有日期过滤器的两步模式
1. 搜索机构获取ID
2. 按机构ID和年份过滤作品

**Python 示例**：
<<<代码块_1>>>

## 某一主题的高被引论文

**用户查询**：“查找过去 5 年有关 CRISPR 被引用最多的论文”

**方法**：搜索+过滤+排序

**Python 示例**：
<<<代码块_2>>>

## 关于某个主题的开放获取论文

**用户查询**：“查找有关气候变化的开放获取论文”

**方法**：搜索+OA过滤

**Python 示例**：
<<<代码块_3>>>

## 出版趋势分析

**用户查询**：“显示多年来机器学习的出版趋势”

**方法**：使用group_by按年份聚合

**Python 示例**：
<<<代码块_4>>>

## 分析研究成果

**用户查询**：“斯坦福大学2020-2024年研究成果分析”

**方法**：多重聚合进行综合分析

**Python 示例**：
<<<代码块_5>>>

## 按 DOI 查找论文（批量）

**用户查询**：“获取这 10 个 DOI 的信息：...”

**方法**：使用管道分隔符批量查找

**Python 示例**：
<<<代码块_6>>>

## 论文随机样本

**用户查询**：“给我 2023 年的 50 篇随机论文”

**方法**：使用带有种子的样本参数来实现可重复性

**Python 示例**：
```python
works = client.sample_works(
    sample_size=50,
    seed=42,  # For reproducibility
    filter_params={
        "publication_year": "2023",
        "is_oa": "true"
    }
)

print(f"Got {len(works)} random papers from 2023")
```

## 多个机构的论文

**用户查询**：“查找麻省理工学院和斯坦福大学作者的论文”

**方法**：在同一属性内使用 + 运算符进行 AND

**Python 示例**：
```python
# First, get institution IDs
mit_response = client._make_request(
    '/institutions',
    params={'search': 'MIT', 'per-page': 1}
)
mit_id = mit_response['results'][0]['id'].split('/')[-1]

stanford_response = client._make_request(
    '/institutions',
    params={'search': 'Stanford', 'per-page': 1}
)
stanford_id = stanford_response['results'][0]['id'].split('/')[-1]

# Find works with authors from both institutions
works = client.search_works(
    filter_params={
        "authorships.institutions.id": f"{mit_id}+{stanford_id}"
    },
    per_page=100
)

print(f"Found {works['meta']['count']} collaborative papers")
```

## 特定期刊中的论文

**用户查询**：“获取《Nature》2023年发表的所有论文”

**方法**：两步 - 找到期刊ID，然后过滤作品

**Python 示例**：
```python
# Step 1: Find journal source ID
source_response = client._make_request(
    '/sources',
    params={'search': 'Nature', 'per-page': 1}
)
source = source_response['results'][0]
source_id = source['id'].split('/')[-1]

print(f"Found journal: {source['display_name']} (ID: {source_id})")

# Step 2: Get works from that source
works = client.search_works(
    filter_params={
        "primary_location.source.id": source_id,
        "publication_year": "2023"
    },
    per_page=200
)

print(f"Found {works['meta']['count']} papers from Nature in 2023")
```

## 按机构划分的主题分析

**用户查询**：“麻省理工学院研究最多的主题是什么？”

**方法**：按机构过滤，按主题分组

**Python 示例**：
```python
# Get MIT ID
inst_response = client._make_request(
    '/institutions',
    params={'search': 'MIT', 'per-page': 1}
)
mit_id = inst_response['results'][0]['id'].split('/')[-1]

# Group by topics
topics = client.group_by(
    entity_type='works',
    group_field='topics.id',
    filter_params={
        "authorships.institutions.id": mit_id,
        "publication_year": ">2020"
    }
)

print("Top research topics at MIT (2020+):")
for i, topic in enumerate(topics[:10], 1):
    print(f"{i}. {topic['key_display_name']}: {topic['count']} works")
```

## 引文分析

**用户查询**：“查找引用此特定 DOI 的论文”

**方法**：通过 DOI 获取工作，然后使用cited_by_api_url

**Python 示例**：
```python
# Get the work
doi = "https://doi.org/10.1038/s41586-021-03819-2"
work = client.get_entity('works', doi)

# Get papers that cite it
cited_by_url = work['cited_by_api_url']

# Extract just the query part and use it
import requests
response = requests.get(cited_by_url, params={'mailto': client.email})
citing_works = response.json()

print(f"{work['title']}")
print(f"Total citations: {work['cited_by_count']}")
print(f"\nRecent citing papers:")
for citing_work in citing_works['results'][:5]:
    print(f"  - {citing_work['title']} ({citing_work['publication_year']})")
```

## 大规模数据提取

**用户查询**：“获取最近3年所有有关量子计算的论文”

**方法**：对所有结果进行分页

**Python 示例**：
```python
all_papers = client.paginate_all(
    endpoint='/works',
    params={
        'search': 'quantum computing',
        'filter': 'publication_year:2022-2024'
    },
    max_results=10000  # Limit to prevent excessive API calls
)

print(f"Retrieved {len(all_papers)} papers")

# Save to CSV
import csv
with open('quantum_papers.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Title', 'Year', 'Citations', 'DOI', 'OA Status'])

    for paper in all_papers:
        writer.writerow([
            paper['title'],
            paper['publication_year'],
            paper['cited_by_count'],
            paper.get('doi', 'N/A'),
            paper['open_access']['oa_status']
        ])
```

## 复杂的多过滤查询

**用户查询**：“查找来自顶级机构的最新、高引用、开放获取的人工智能论文”

**方法**：组合多个过滤器

**Python 示例**：
```python
# Get IDs for top institutions
top_institutions = ['MIT', 'Stanford', 'Oxford']
inst_ids = []

for inst_name in top_institutions:
    response = client._make_request(
        '/institutions',
        params={'search': inst_name, 'per-page': 1}
    )
    if response['results']:
        inst_id = response['results'][0]['id'].split('/')[-1]
        inst_ids.append(inst_id)

# Combine with pipe for OR
inst_filter = '|'.join(inst_ids)

# Complex query
works = client.search_works(
    search="artificial intelligence",
    filter_params={
        "publication_year": ">2022",
        "cited_by_count": ">50",
        "is_oa": "true",
        "authorships.institutions.id": inst_filter
    },
    sort="cited_by_count:desc",
    per_page=200
)

print(f"Found {works['meta']['count']} papers matching criteria")
for work in works['results'][:10]:
    print(f"{work['title']}")
    print(f"  Citations: {work['cited_by_count']}, Year: {work['publication_year']}")
```