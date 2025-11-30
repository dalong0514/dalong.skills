<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： 临床试验数据库
描述：“通过 API v2 查询 ClinicalTrials.gov。按条件、药物、位置、状态或阶段搜索试验。通过 NCT ID 检索试验详细信息，导出数据，用于临床研究和患者匹配。”
---

# ClinicalTrials.gov 数据库

## 概述

ClinicalTrials.gov 是全球临床研究的综合注册中心，由美国国家医学图书馆维护。访问 API v2 来搜索试验、检索详细的研究信息、按各种标准进行过滤以及导出数据进行分析。该 API 是公开的（无需身份验证），速率限制为每分钟约 50 个请求，支持 JSON 和 CSV 格式。

## 何时使用此技能

在以下场景中处理临床试验数据时应使用此技能：

- **患者匹配** - 寻找针对特定情况或患者群体的招募试验
- **研究分析** - 分析临床试验趋势、结果或研究设计
- **药物/干预研究** - 确定测试特定药物或干预措施的试验
- **地理搜索** - 在特定地点或地区定位试验
- **赞助商/组织跟踪** - 查找特定机构进行的试验
- **数据导出** - 提取临床试验数据以进行进一步分析或报告
- **试验监控** - 跟踪特定试验的状态更新或结果
- **资格筛选** - 审查试验的纳入/排除标准

## 快速入门

### 基本搜索查询

使用帮助脚本搜索临床试验：

```bash
cd scientific-databases/clinicaltrials-database/scripts
python3 query_clinicaltrials.py
```

或者直接将 Python 与 `requests` 库一起使用：

<<<代码块_1>>>

### 检索特定试验

使用 NCT ID 获取有关试验的详细信息：

<<<代码块_2>>>

## 核心能力

### 1. 按病情/疾病搜索

使用 `query.cond` 参数查找研究特定医疗状况或疾病的试验。

**示例：查找招募糖尿病试验**

<<<代码块_3>>>

**常见用例：**
- 寻找罕见疾病的试验
- 确定合并症的试验
- 跟踪特定诊断的试验可用性

### 2. 按干预/药物搜索

使用 `query.intr` 参数搜索测试特定干预措施、药物、设备或程序的试验。

**示例：查找测试 Pembrolizumab 的 3 期试验**

<<<代码块_4>>>

**常见用例：**
- 药物开发跟踪
- 制药公司的竞争情报
- 临床医生的治疗选择研究

### 3. 地理搜索

使用 `query.locn` 参数查找特定位置的试验。

**示例：查找纽约的癌症试验**

<<<代码块_5>>>

**常见用例：**
- 患者转诊至当地试验
- 试验地理分布分析
- 新试验的地点选择

### 4. 按赞助商/组织搜索

使用 `query.spons` 参数查找特定组织进行的试验。

**示例：查找 NCI 赞助的试验**

<<<代码块_6>>>

**常见用例：**
- 跟踪机构研究组合
- 分析资助组织的优先事项
- 确定合作机会

### 5. 按研究状态过滤

使用 `filter.overallStatus` 参数按招募或完成状态筛选试验。

**有效状态值：**
- `RECRUITING` - 目前正在招募参与者
- `NOT_YET_RECRUITING` - 尚未开放招募
- `ENROLLING_BY_INVITATION` - 仅通过邀请注册
- `ACTIVE_NOT_RECRUITING` - 活跃但不再招募
- `SUSPENDED` - 暂时停止
- `TERMINATED` - 过早停止
- `COMPLETED` - 研究已结束
- `WITHDRAWN` - 注册前撤回

**示例：查找最近完成的试验及结果**

```python
from scripts.query_clinicaltrials import search_studies

results = search_studies(
    condition="alzheimer disease",
    status="COMPLETED",
    sort="LastUpdatePostDate:desc",
    page_size=50
)

# Filter for trials with results
trials_with_results = [
    study for study in results['studies']
    if study.get('hasResults', False)
]

print(f"Found {len(trials_with_results)} completed trials with results")
```

### 6. 检索详细的研究信息

获取有关特定试验的全面信息，包括资格标准、结果、联系人和地点。

**示例：提取资格标准**

```python
from scripts.query_clinicaltrials import get_study_details

study = get_study_details("NCT04852770")
eligibility = study['protocolSection']['eligibilityModule']

print(f"Eligible Ages: {eligibility.get('minimumAge')} - {eligibility.get('maximumAge')}")
print(f"Eligible Sex: {eligibility.get('sex')}")
print(f"\nInclusion Criteria:")
print(eligibility.get('eligibilityCriteria'))
```

**示例：提取联系信息**

```python
from scripts.query_clinicaltrials import get_study_details

study = get_study_details("NCT04852770")
contacts_module = study['protocolSection']['contactsLocationsModule']

# Overall contacts
if 'centralContacts' in contacts_module:
    for contact in contacts_module['centralContacts']:
        print(f"Contact: {contact.get('name')}")
        print(f"Phone: {contact.get('phone')}")
        print(f"Email: {contact.get('email')}")

# Study locations
if 'locations' in contacts_module:
    for location in contacts_module['locations']:
        print(f"\nFacility: {location.get('facility')}")
        print(f"City: {location.get('city')}, {location.get('state')}")
        if location.get('status'):
            print(f"Status: {location['status']}")
```

### 7. 分页和批量数据检索

使用分页有效处理大型结果集。

**示例：检索所有匹配的试验**

```python
from scripts.query_clinicaltrials import search_with_all_results

# Get all trials (automatically handles pagination)
all_trials = search_with_all_results(
    condition="rare disease",
    status="RECRUITING"
)

print(f"Retrieved {len(all_trials)} total trials")
```

**示例：带控件的手动分页**

```python
from scripts.query_clinicaltrials import search_studies

all_studies = []
page_token = None
max_pages = 10  # Limit to avoid excessive requests

for page in range(max_pages):
    results = search_studies(
        condition="cancer",
        page_size=1000,  # Max page size
        page_token=page_token
    )

    all_studies.extend(results['studies'])

    # Check for next page
    page_token = results.get('pageToken')
    if not page_token:
        break

print(f"Retrieved {len(all_studies)} studies across {page + 1} pages")
```

### 8. 数据导出至 CSV

将试验数据导出为 CSV 格式，以便在电子表格软件或数据分析工具中进行分析。
**示例：导出到 CSV 文件**

```python
from scripts.query_clinicaltrials import search_studies

# Request CSV format
results = search_studies(
    condition="heart disease",
    status="RECRUITING",
    format="csv",
    page_size=1000
)

# Save to file
with open("heart_disease_trials.csv", "w") as f:
    f.write(results)

print("Data exported to heart_disease_trials.csv")
```

**注意：** CSV 格式返回字符串而不是 JSON 字典。

### 9. 提取和总结研究信息

提取关键信息以进行快速概述或报告。

**示例：创建试用摘要**

```python
from scripts.query_clinicaltrials import get_study_details, extract_study_summary

# Get details and extract summary
study = get_study_details("NCT04852770")
summary = extract_study_summary(study)

print(f"NCT ID: {summary['nct_id']}")
print(f"Title: {summary['title']}")
print(f"Status: {summary['status']}")
print(f"Phase: {', '.join(summary['phase'])}")
print(f"Enrollment: {summary['enrollment']}")
print(f"Last Update: {summary['last_update']}")
print(f"\nBrief Summary:\n{summary['brief_summary']}")
```

### 10.组合查询策略

组合多个过滤器以进行有针对性的搜索。

**示例：多标准搜索**

```python
from scripts.query_clinicaltrials import search_studies

# Find Phase 2/3 immunotherapy trials for lung cancer in California
results = search_studies(
    condition="lung cancer",
    intervention="immunotherapy",
    location="California",
    status=["RECRUITING", "NOT_YET_RECRUITING"],
    page_size=100
)

# Further filter by phase
phase2_3_trials = [
    study for study in results['studies']
    if any(phase in ['PHASE2', 'PHASE3']
           for phase in study['protocolSection'].get('designModule', {}).get('phases', []))
]

print(f"Found {len(phase2_3_trials)} Phase 2/3 immunotherapy trials")
```

## 资源

### 脚本/query_clinicaltrials.py

综合 Python 脚本为常见查询模式提供辅助函数：

- `search_studies()` - 使用各种过滤器搜索试验
- `get_study_details()` - 检索特定试验的完整信息
- `search_with_all_results()` - 自动对所有结果进行分页
- `extract_study_summary()` - 提取关键信息以快速概览

直接运行脚本，示例用法：

```bash
python3 scripts/query_clinicaltrials.py
```

### 参考文献/api_reference.md

详细的 API 文档包括：

- 完整的端点规格
- 所有查询参数和有效值
- 响应数据结构和模块
- 常见用例和代码示例
- 错误处理和最佳实践
- 数据标准（ISO 8601 日期、CommonMark markdown）

使用不熟悉的 API 功能或排除问题时加载此参考。

## 最佳实践

### 速率限制管理

API 的速率限制约为每分钟 50 个请求。对于批量数据检索：

1.使用最大页面大小（1000）来最小化请求
2. 对速率限制错误（429状态）实施指数退避
3. 在大规模数据收集的请求之间添加延迟

```python
import time
import requests

def search_with_rate_limit(params):
    try:
        response = requests.get("https://clinicaltrials.gov/api/v2/studies", params=params)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 429:
            print("Rate limited. Waiting 60 seconds...")
            time.sleep(60)
            return search_with_rate_limit(params)  # Retry
        raise
```

### 数据结构导航

API 响应具有嵌套结构。通用信息的关键路径：

- **NCT ID**：`study['protocolSection']['identificationModule']['nctId']`
- **标题**：`study['protocolSection']['identificationModule']['briefTitle']`
- **状态**：`study['protocolSection']['statusModule']['overallStatus']`
- **阶段**：`study['protocolSection']['designModule']['phases']`
- **资格**：`study['protocolSection']['eligibilityModule']`
- **位置**：`study['protocolSection']['contactsLocationsModule']['locations']`
- **干预**：`study['protocolSection']['armsInterventionsModule']['interventions']`

### 错误处理

始终对网络请求实施正确的错误处理：

```python
import requests

try:
    response = requests.get(url, params=params, timeout=30)
    response.raise_for_status()
    data = response.json()
except requests.exceptions.HTTPError as e:
    print(f"HTTP error: {e.response.status_code}")
except requests.exceptions.RequestException as e:
    print(f"Request failed: {e}")
except ValueError as e:
    print(f"JSON decode error: {e}")
```

### 处理缺失数据

并非所有试验都有完整的信息。始终检查字段是否存在：

```python
# Safe navigation with .get()
phases = study['protocolSection'].get('designModule', {}).get('phases', [])
enrollment = study['protocolSection'].get('designModule', {}).get('enrollmentInfo', {}).get('count', 'N/A')

# Check before accessing
if 'resultsSection' in study:
    # Process results
    pass
```

## 技术规格

- **基本 URL**：`https://clinicaltrials.gov/api/v2`
- **身份验证**：不需要（公共 API）
- **速率限制**：每个 IP 约 50 个请求/分钟
- **响应格式**：JSON（默认）、CSV
- **最大页面大小**：每个请求 1000 个研究
- **日期格式**：ISO 8601
- **文本格式**：用于富文本字段的 CommonMark Markdown
- **API 版本**：2.0（2024 年 3 月发布）
- **API 规范**：OpenAPI 3.0

有关完整的技术详细信息，请参阅`references/api_reference.md`。