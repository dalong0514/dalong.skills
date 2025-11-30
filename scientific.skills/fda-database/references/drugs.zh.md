<!-- 此文件由机器翻译自 drugs.md -->

# FDA 药物数据库

该参考文献涵盖了可通过 openFDA 访问的所有 FDA 药物相关 API 端点。

## 概述

FDA 药物数据库提供有关药品的信息，包括不良事件、标签、召回、批准和短缺。所有端点都遵循 openFDA API 结构并返回 JSON 格式的数据。

## 可用端点

### 1. 药物不良事件

**端点**：`https://api.fda.gov/drug/event.json`

**目的**：获取提交给 FDA 的药物副作用、产品使用错误、产品质量问题和治疗失败的报告。

**数据来源**：FDA 不良事件报告系统 (FAERS)

**关键字段**：
- `patient.drug.medicinalproduct` - 药物名称
- `patient.drug.drugindication` - 服药原因
- `patient.reaction.reactionmeddrapt` - 不良反应描述
- `receivedate` - 收到报告的日期
- `serious` - 事件是否严重（1 = 严重，2 = 不严重）
- `seriousnessdeath` - 该事件是否导致死亡
- `primarysource.qualification` - 报告者资格（医生、药剂师等）

**常见用例**：
- 安全信号检测
- 上市后监督
- 药物相互作用分析
- 比较安全性研究

**查询示例**：
```python
# Find adverse events for a specific drug
import requests

api_key = "YOUR_API_KEY"
url = "https://api.fda.gov/drug/event.json"

# Search for aspirin-related adverse events
params = {
    "api_key": api_key,
    "search": "patient.drug.medicinalproduct:aspirin",
    "limit": 10
}

response = requests.get(url, params=params)
data = response.json()
```

<<<代码块_1>>>

### 2. 药品标签

**端点**：`https://api.fda.gov/drug/label.json`

**目的**：访问结构化产品信息，包括 FDA 批准和销售的药品的处方信息、警告、适应症和用途。

**数据来源**：结构化产品标签 (SPL)

**关键字段**：
- `openfda.brand_name` - 药物的品牌名称
- `openfda.generic_name` - 通用名称
- `indications_and_usage` - 批准的用途
- `warnings` - 重要安全警告
- `adverse_reactions` - 已知不良反应
- `dosage_and_administration` - 如何使用药物
- `description` - 化学和物理描述
- `pharmacodynamics` - 药物的作用原理
- `contraindications` - 何时不使用该药物
- `drug_interactions` - 已知药物相互作用
- `active_ingredient` - 活性成分
- `inactive_ingredient` - 非活性成分

**常见用例**：
- 临床决策支持
- 药品信息查询
- 患者教育材料
- 处方管理
- 药物比较分析

**查询示例**：
<<<代码块_2>>>

<<<代码块_3>>>

### 3. 国家药品法规 (NDC) 目录

**端点**：`https://api.fda.gov/drug/ndc.json`

**目的**：访问 NDC 目录，其中包含有关国家药品代码识别的药品信息。

**数据来源**：FDA NDC 目录

**关键字段**：
- `product_ndc` - 10 位 NDC 产品标识符
- `generic_name` - 通用药物名称
- `labeler_name` - 制造/分销的公司
- `brand_name` - 品牌名称（如果适用）
- `dosage_form` - 剂型（片剂、胶囊、溶液等）
- `route` - 给药途径（口服、注射、外用等）
- `product_type` - 药品类型
- `marketing_category` - 监管途径（NDA、ANDA、OTC 等）
- `application_number` - FDA 申请号
- `active_ingredients` - 具有优势的活性成分列表
- `packaging` - 包描述和 NDC 代码
- `listing_expiration_date` - 列表过期时

**常见用例**：
- NDC查找和验证
- 产品标识
- 供应链管理
- 处方处理
- 保险索赔处理

**查询示例**：
<<<代码块_4>>>

<<<代码块_5>>>

<<<代码块_6>>>

### 4. 药品召回执行报告

**端点**：`https://api.fda.gov/drug/enforcement.json`

**目的**：获取 FDA 发布的药品召回执行报告。

**数据来源**：FDA 执法报告

**关键字段**：
- `status` - 当前状态（正在进行、已完成、已终止）
- `recall_number` - 唯一召回标识符
- `classification` - I、II 或 III 级（严重性）
- `product_description` - 召回产品的描述
- `reason_for_recall` - 产品为何被召回
- `product_quantity` - 召回的产品数量
- `code_info` - 批号、序列号、NDC
- `distribution_pattern` - 地理分布
- `recalling_firm` - 公司进行召回
- `recall_initiation_date` - 召回开始时
- `report_date` - 当 FDA 收到通知时
- `voluntary_mandated` - 召回类型

**分类级别**：
- **I 类**：可能导致严重健康问题或死亡的危险或有缺陷的产品
- **II 类**：可能导致暂时健康问题或构成轻微严重威胁的产品
- **III 类**：不太可能引起不良健康反应但违反 FDA 标签/制造规定的产品

**常见用例**：
- 质量保证监控
- 供应链风险管理
- 患者安全警报
- 监管合规性跟踪

**查询示例**：
```python
# Find all Class I (most serious) drug recalls
params = {
    "api_key": api_key,
    "search": "classification:Class+I",
    "limit": 20,
    "sort": "report_date:desc"
}

response = requests.get("https://api.fda.gov/drug/enforcement.json", params=params)
```

```python
# Search for recalls of a specific drug
params = {
    "api_key": api_key,
    "search": "product_description:*metformin*",
    "limit": 10
}
```

```python
# Find ongoing recalls
params = {
    "api_key": api_key,
    "search": "status:Ongoing",
    "limit": 50
}
```

### 5. 药物@FDA

**端点**：`https://api.fda.gov/drug/drugsfda.json`

**目的**：从 Drugs@FDA 数据库获取 FDA 批准的药品的全面信息，包括批准历史和监管信息。

**数据来源**：Drugs@FDA 数据库（自 1939 年以来批准的大多数药物）

**关键字段**：
- `application_number` - NDA/ANDA/BLA 编号
- `sponsor_name` - 提交申请的公司
- `openfda.brand_name` - 品牌名称
- `openfda.generic_name` - 通用名称
- `products` - 本申请下批准的产品数组
- `products.active_ingredients` - 具有优势的活性成分
- `products.dosage_form` - 剂型
- `products.route` - 给药途径
- `products.marketing_status` - 当前营销状态
- `submissions` - 监管提交数组
- `submissions.submission_type` - 提交类型
- `submissions.submission_status` - 状态（已批准、待处理等）
- `submissions.submission_status_date` - 状态日期
- `submissions.review_priority` - 优先或标准审核

**常见用例**：
- 药品审批研究
- 监管途径分析
- 历史审批跟踪
- 竞争情报
- 市场准入研究

**查询示例**：
```python
# Find approval information for a specific drug
params = {
    "api_key": api_key,
    "search": "openfda.brand_name:Keytruda",
    "limit": 1
}

response = requests.get("https://api.fda.gov/drug/drugsfda.json", params=params)
```

```python
# Get all drugs approved by a specific sponsor
params = {
    "api_key": api_key,
    "search": "sponsor_name:Moderna",
    "limit": 100
}
```

```python
# Find drugs with priority review designation
params = {
    "api_key": api_key,
    "search": "submissions.review_priority:Priority",
    "limit": 50
}
```

### 6. 药品短缺

**端点**：`https://api.fda.gov/drug/drugshortages.json`

**目的**：获取有关影响美国的当前和已解决的药物短缺的信息。

**数据来源**：FDA药品短缺数据库

**关键字段**：
- `product_name` - 短缺药品名称
- `status` - 当前状态（当前短缺、已解决、已停产）
- `reason` - 短缺原因
- `shortage_start_date` - 短缺开始时
- `resolution_date` - 短缺问题得到解决时（如果适用）
- `discontinuation_date` - 如果产品已停产
- `active_ingredient` - 活性成分
- `marketed_by` - 营销该产品的公司
- `presentation` - 剂型和强度

**常见用例**：
- 处方管理
- 供应链规划
- 病人护理的连续性
- 治疗替代鉴定
- 采购计划

**查询示例**：
```python
# Find current drug shortages
params = {
    "api_key": api_key,
    "search": "status:Currently+in+Shortage",
    "limit": 100
}

response = requests.get("https://api.fda.gov/drug/drugshortages.json", params=params)
```

```python
# Search for shortages of a specific drug
params = {
    "api_key": api_key,
    "search": "product_name:*amoxicillin*",
    "limit": 10
}
```

```python
# Get shortage history (both current and resolved)
params = {
    "api_key": api_key,
    "search": "active_ingredient:epinephrine",
    "limit": 50
}
```

## 集成技巧

### 错误处理

```python
import requests
import time

def query_fda_drug(endpoint, params, max_retries=3):
    """
    Query FDA drug database with error handling and retry logic.

    Args:
        endpoint: Full URL endpoint (e.g., "https://api.fda.gov/drug/event.json")
        params: Dictionary of query parameters
        max_retries: Maximum number of retry attempts

    Returns:
        Response JSON data or None if error
    """
    for attempt in range(max_retries):
        try:
            response = requests.get(endpoint, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError as e:
            if response.status_code == 404:
                print(f"No results found for query")
                return None
            elif response.status_code == 429:
                # Rate limit exceeded, wait and retry
                wait_time = 60 * (attempt + 1)
                print(f"Rate limit exceeded. Waiting {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print(f"HTTP error occurred: {e}")
                return None
        except requests.exceptions.RequestException as e:
            print(f"Request error: {e}")
            if attempt < max_retries - 1:
                time.sleep(5)
            else:
                return None
    return None
```

### 大型结果集的分页

```python
def get_all_results(endpoint, search_query, api_key, max_results=1000):
    """
    Retrieve all results for a query using pagination.

    Args:
        endpoint: API endpoint URL
        search_query: Search query string
        api_key: FDA API key
        max_results: Maximum total results to retrieve

    Returns:
        List of all result records
    """
    all_results = []
    skip = 0
    limit = 100  # Max per request

    while len(all_results) < max_results:
        params = {
            "api_key": api_key,
            "search": search_query,
            "limit": limit,
            "skip": skip
        }

        data = query_fda_drug(endpoint, params)
        if not data or "results" not in data:
            break

        results = data["results"]
        all_results.extend(results)

        # Check if we've retrieved all available results
        if len(results) < limit:
            break

        skip += limit
        time.sleep(0.25)  # Rate limiting courtesy

    return all_results[:max_results]
```

## 最佳实践

1. **始终使用 HTTPS** - 不接受 HTTP 请求
2. **包括 API 密钥** - 提供更高的速率限制（120,000/天 vs 1,000/天）
3. **使用精确匹配进行聚合** - 在计数查询中为字段名称添加 `.exact` 后缀
4. **实施速率限制** - 保持在 240 个请求/分钟以内
5. **缓存结果** - 避免对相同数据进行冗余查询
6. **优雅地处理错误** - 实现瞬态故障的重试逻辑
7. **使用特定字段搜索** - 比全文搜索更高效
8. **验证 NDC 代码** - 使用标准 11 位数字格式并删除连字符
9. **监控 API 状态** - 检查 openFDA 状态页面是否出现中断
10. **尊重数据限制** - OpenFDA 仅包含公共数据，并非所有 FDA 数据

## 其他资源

- OpenFDA 药物 API 文档：https://open.fda.gov/apis/drug/
- API 基础知识：请参阅此参考目录中的 `api_basics.md`
- Python 示例：参见 `scripts/fda_drug_query.py`
- 现场参考指南：可在 https://open.fda.gov/apis/drug/[endpoint]/searchable-fields/ 获取