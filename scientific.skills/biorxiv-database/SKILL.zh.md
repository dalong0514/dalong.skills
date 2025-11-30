<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：biorxiv-数据库
描述：bioRxiv 预印本服务器的高效数据库搜索工具。在按关键字、作者、日期范围或类别搜索生命科学预印本、检索论文元数据、下载 PDF 或进行文献综述时，请使用此技能。
---

#bioRxiv 数据库

## 概述

这项技能提供了基于 Python 的高效工具，用于从 bioRxiv 数据库中搜索和检索预印本。它支持按关键字、作者、日期范围和类别进行全面搜索，返回结构化 JSON 元数据，包括标题、摘要、DOI 和引文信息。该技能还支持 PDF 下载以进行全文分析。

## 何时使用此技能

在以下情况下使用此技能：
- 搜索特定研究领域的最新预印本
- 跟踪特定作者的出版物
- 进行系统的文献综述
- 分析一段时间内的研究趋势
- 检索元数据以进行引文管理
- 下载预印本 PDF 进行分析
- 按bioRxiv主题类别过滤论文

## 核心搜索功能

### 1. 关键字搜索

搜索标题、摘要或作者列表中包含特定关键字的预印本。

**基本用法：**
```python
python scripts/biorxiv_search.py \
  --keywords "CRISPR" "gene editing" \
  --start-date 2024-01-01 \
  --end-date 2024-12-31 \
  --output results.json
```

**使用类别过滤器：**
<<<代码块_1>>>

**搜索字段：**
默认情况下，会在标题和摘要中搜索关键字。使用 `--search-fields` 进行自定义：
<<<代码块_2>>>

### 2. 作者搜索

查找特定作者在某个日期范围内的所有论文。

**基本用法：**
<<<代码块_3>>>

**最近发表的文章：**
<<<代码块_4>>>

### 3. 日期范围搜索

检索特定日期范围内发布的所有预印本。

**基本用法：**
<<<代码块_5>>>

**使用类别过滤器：**
<<<代码块_6>>>

**回溯快捷方式：**
```python
# Last 30 days
python scripts/biorxiv_search.py \
  --days-back 30 \
  --output last_month.json
```

### 4. DOI 的论文详细信息

检索特定预印本的详细元数据。

**基本用法：**
```python
python scripts/biorxiv_search.py \
  --doi "10.1101/2024.01.15.123456" \
  --output paper_details.json
```

**接受完整的 DOI URL：**
```python
python scripts/biorxiv_search.py \
  --doi "https://doi.org/10.1101/2024.01.15.123456"
```

### 5.PDF下载

下载任何预印本的全文 PDF。

**基本用法：**
```python
python scripts/biorxiv_search.py \
  --doi "10.1101/2024.01.15.123456" \
  --download-pdf paper.pdf
```

**批量处理：**
对于多个 PDF，从搜索结果 JSON 中提取 DOI 并下载每篇论文：
```python
import json
from biorxiv_search import BioRxivSearcher

# Load search results
with open('results.json') as f:
    data = json.load(f)

searcher = BioRxivSearcher(verbose=True)

# Download each paper
for i, paper in enumerate(data['results'][:10]):  # First 10 papers
    doi = paper['doi']
    searcher.download_pdf(doi, f"papers/paper_{i+1}.pdf")
```

## 有效类别

按 bioRxiv 主题类别过滤搜索：

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

## 输出格式

所有搜索都会返回具有以下格式的结构化 JSON：

```json
{
  "query": {
    "keywords": ["CRISPR"],
    "start_date": "2024-01-01",
    "end_date": "2024-12-31",
    "category": "genomics"
  },
  "result_count": 42,
  "results": [
    {
      "doi": "10.1101/2024.01.15.123456",
      "title": "Paper Title Here",
      "authors": "Smith J, Doe J, Johnson A",
      "author_corresponding": "Smith J",
      "author_corresponding_institution": "University Example",
      "date": "2024-01-15",
      "version": "1",
      "type": "new results",
      "license": "cc_by",
      "category": "genomics",
      "abstract": "Full abstract text...",
      "pdf_url": "https://www.biorxiv.org/content/10.1101/2024.01.15.123456v1.full.pdf",
      "html_url": "https://www.biorxiv.org/content/10.1101/2024.01.15.123456v1",
      "jatsxml": "https://www.biorxiv.org/content/...",
      "published": ""
    }
  ]
}
```

## 常见使用模式

### 文献综述工作流程

1. **广泛的关键词搜索：**
```python
python scripts/biorxiv_search.py \
  --keywords "organoids" "tissue engineering" \
  --start-date 2023-01-01 \
  --end-date 2024-12-31 \
  --category bioengineering \
  --output organoid_papers.json
```

2. **提取并查看结果：**
```python
import json

with open('organoid_papers.json') as f:
    data = json.load(f)

print(f"Found {data['result_count']} papers")

for paper in data['results'][:5]:
    print(f"\nTitle: {paper['title']}")
    print(f"Authors: {paper['authors']}")
    print(f"Date: {paper['date']}")
    print(f"DOI: {paper['doi']}")
```

3. **下载所选论文：**
```python
from biorxiv_search import BioRxivSearcher

searcher = BioRxivSearcher()
selected_dois = ["10.1101/2024.01.15.123456", "10.1101/2024.02.20.789012"]

for doi in selected_dois:
    filename = doi.replace("/", "_").replace(".", "_") + ".pdf"
    searcher.download_pdf(doi, f"papers/{filename}")
```

### 趋势分析

通过分析一段时间内的发表频率来跟踪研究趋势：

```python
python scripts/biorxiv_search.py \
  --keywords "machine learning" \
  --start-date 2020-01-01 \
  --end-date 2024-12-31 \
  --category bioinformatics \
  --output ml_trends.json
```

然后分析结果的时间分布。

### 作者追踪

监控特定研究人员的预印本：

```python
# Track multiple authors
authors = ["Smith", "Johnson", "Williams"]

for author in authors:
    python scripts/biorxiv_search.py \
      --author "{author}" \
      --days-back 365 \
      --output "{author}_papers.json"
```

## Python API 使用

对于更复杂的工作流程，直接导入并使用 `BioRxivSearcher` 类：

```python
from scripts.biorxiv_search import BioRxivSearcher

# Initialize
searcher = BioRxivSearcher(verbose=True)

# Multiple search operations
keywords_papers = searcher.search_by_keywords(
    keywords=["CRISPR", "gene editing"],
    start_date="2024-01-01",
    end_date="2024-12-31",
    category="genomics"
)

author_papers = searcher.search_by_author(
    author_name="Smith",
    start_date="2023-01-01",
    end_date="2024-12-31"
)

# Get specific paper details
paper = searcher.get_paper_details("10.1101/2024.01.15.123456")

# Download PDF
success = searcher.download_pdf(
    doi="10.1101/2024.01.15.123456",
    output_path="paper.pdf"
)

# Format results consistently
formatted = searcher.format_result(paper, include_abstract=True)
```

## 最佳实践

1. **使用适当的日期范围**：较小的日期范围返回速度更快。对于长时间的关键字搜索，请考虑拆分为多个查询。

2. **按类别过滤**：如果可能，使用`--category`来减少数据传输并提高搜索精度。

3. **尊重率限制**：脚本包含自动延迟（请求之间 0.5 秒）。对于大规模数据收集，请添加额外的延迟。

4. **缓存结果**：将搜索结果保存到JSON文件中，以避免重复的API调用。

5. **版本跟踪**：预印本可以有多个版本。 `version` 字段指示返回哪个版本。 PDF URL 包含版本号。

6. **优雅地处理错误**：检查输出 JSON 中的 `result_count`。空结果可能表示日期范围问题或 API 连接问题。
7. **用于调试的详细模式**：使用 `--verbose` 标志查看 API 请求和响应的详细日志记录。

## 高级功能

### 自定义日期范围逻辑

```python
from datetime import datetime, timedelta

# Last quarter
end_date = datetime.now()
start_date = end_date - timedelta(days=90)

python scripts/biorxiv_search.py \
  --start-date {start_date.strftime('%Y-%m-%d')} \
  --end-date {end_date.strftime('%Y-%m-%d')}
```

### 结果限制

限制返回结果的数量：

```python
python scripts/biorxiv_search.py \
  --keywords "COVID-19" \
  --days-back 30 \
  --limit 50 \
  --output covid_top50.json
```

### 为了速度而排除摘要

当仅需要元数据时：

```python
# Note: Abstract inclusion is controlled in Python API
from scripts.biorxiv_search import BioRxivSearcher

searcher = BioRxivSearcher()
papers = searcher.search_by_keywords(keywords=["AI"], days_back=30)
formatted = [searcher.format_result(p, include_abstract=False) for p in papers]
```

## 程序化集成

将搜索结果集成到下游分析管道中：

```python
import json
import pandas as pd

# Load results
with open('results.json') as f:
    data = json.load(f)

# Convert to DataFrame for analysis
df = pd.DataFrame(data['results'])

# Analyze
print(f"Total papers: {len(df)}")
print(f"Date range: {df['date'].min()} to {df['date'].max()}")
print(f"\nTop authors by paper count:")
print(df['authors'].str.split(',').explode().str.strip().value_counts().head(10))

# Filter and export
recent = df[df['date'] >= '2024-06-01']
recent.to_csv('recent_papers.csv', index=False)
```

## 测试技能

要验证 bioRxiv 数据库技能是否正常工作，请运行综合测试套件。

**先决条件：**
```bash
uv pip install requests
```

**运行测试：**
```bash
python tests/test_biorxiv_search.py
```

测试套件验证：
- **初始化**：BioRxivSearcher 类实例化
- **日期范围搜索**：检索特定日期范围内的论文
- **类别过滤**：按bioRxiv类别过滤论文
- **关键字搜索**：查找包含特定关键字的论文
- **DOI 查找**：通过 DOI 检索特定论文
- **结果格式**：纸张元数据的正确格式
- **间隔搜索**：按时间间隔获取最近的论文

**预期输出：**
```
🧬 bioRxiv Database Search Skill Test Suite
======================================================================

🧪 Test 1: Initialization
✅ BioRxivSearcher initialized successfully

🧪 Test 2: Date Range Search
✅ Found 150 papers between 2024-01-01 and 2024-01-07
   First paper: Novel CRISPR-based approach for genome editing...

[... additional tests ...]

======================================================================
📊 Test Summary
======================================================================
✅ PASS: Initialization
✅ PASS: Date Range Search
✅ PASS: Category Filtering
✅ PASS: Keyword Search
✅ PASS: DOI Lookup
✅ PASS: Result Formatting
✅ PASS: Interval Search
======================================================================
Results: 7/7 tests passed (100%)
======================================================================

🎉 All tests passed! The bioRxiv database skill is working correctly.
```

**注意：** 如果在特定日期范围或类别中未找到论文，某些测试可能会显示警告。这是正常现象，并不表示出现故障。

## 参考文档

有关详细的 API 规范、端点文档和响应架构，请参阅：
- `references/api_reference.md` - 完整的bioRxiv API文档

参考文件包括：
- 完整的API端点规范
- 响应格式详细信息
- 错误处理模式
- 速率限制指南
- 高级搜索模式