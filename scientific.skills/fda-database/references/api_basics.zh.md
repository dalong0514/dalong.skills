<!-- 此文件由机器翻译自 api_basics.md -->

# OpenFDA API 基础知识

该参考提供了有关使用 openFDA API 的全面信息，包括身份验证、速率限制、查询语法和最佳实践。

## 开始使用

### 基本网址

所有 openFDA API 端点都遵循以下结构：
```
https://api.fda.gov/{category}/{endpoint}.json
```

示例：
- `https://api.fda.gov/drug/event.json`
- `https://api.fda.gov/device/510k.json`
- `https://api.fda.gov/food/enforcement.json`

### 需要 HTTPS

**所有请求必须使用 HTTPS**。 HTTP 请求不被接受并且将会失败。

## 身份验证

### API 密钥注册

虽然 openFDA 可以在没有 API 密钥的情况下使用，但强烈建议注册免费的 API 密钥以获得更高的速率限制。

**注册**：访问https://open.fda.gov/apis/authentication/进行注册

**API 密钥的好处**：
- 更高的速率限制（240 请求/分钟，120,000 请求/天）
- 更适合生产应用
- 无需额外费用

### 使用您的 API 密钥

使用以下两种方法之一将您的 API 密钥包含在请求中：

**方法一：查询参数（推荐）**
<<<代码块_1>>>

**方法2：基本身份验证**
<<<代码块_2>>>

## 速率限制

### 当前限制

|状态 |每分钟请求数 |每日请求数 |
|--------|--------------------|------------------|
| **没有 API 密钥** |每个 IP 地址 240 |每个 IP 地址 1,000 个 |
| **使用 API 密钥** |每把 240 |每把 120,000 个 |

### 速率限制标头

API 在响应标头中返回速率限制信息：
<<<代码块_3>>>

### 处理率限制

当您超过速率限制时，API 将返回：
- **状态代码**：`429 Too Many Requests`
- **错误消息**：表示超出速率限制

**最佳实践**：实施指数退避：
<<<代码块_4>>>

### 增加限制

对于需要更高限制的应用程序，请通过 openFDA 团队的网站联系他们，并提供有关您的用例的详细信息。

## 查询语法

### 基本结构

查询使用以下格式：
<<<代码块_5>>>

参数由 & 符号分隔 (`&`)。

### 搜索参数

`search` 参数是过滤结果的主要方式。

**基本格式**：
<<<代码块_6>>>

**示例**：
```python
params = {
    "api_key": api_key,
    "search": "patient.drug.medicinalproduct:aspirin"
}
```

### 搜索运算符

#### 与运算符
组合多个条件（两者都必须为真）：
```python
# Find aspirin adverse events in Canada
params = {
    "search": "patient.drug.medicinalproduct:aspirin+AND+occurcountry:ca"
}
```

#### 或运算符
任一条件都可以为真（OR 是用空格隐含的）：
```python
# Find aspirin OR ibuprofen
params = {
    "search": "patient.drug.medicinalproduct:(aspirin ibuprofen)"
}
```

或者明确地：
```python
params = {
    "search": "patient.drug.medicinalproduct:aspirin+OR+patient.drug.medicinalproduct:ibuprofen"
}
```

#### NOT 运算符
排除结果：
```python
# Events NOT in the United States
params = {
    "search": "_exists_:occurcountry+AND+NOT+occurcountry:us"
}
```

#### 通配符
使用星号 (`*`) 进行部分匹配：
```python
# Any drug starting with "met"
params = {
    "search": "patient.drug.medicinalproduct:met*"
}

# Any drug containing "cillin"
params = {
    "search": "patient.drug.medicinalproduct:*cillin*"
}
```

#### 精确短语匹配
对确切的短语使用引号：
```python
params = {
    "search": 'patient.reaction.reactionmeddrapt:"heart attack"'
}
```

#### 范围查询
在范围内搜索：
```python
# Date range (YYYYMMDD format)
params = {
    "search": "receivedate:[20200101+TO+20201231]"
}

# Numeric range
params = {
    "search": "patient.patientonsetage:[18+TO+65]"
}

# Open-ended ranges
params = {
    "search": "patient.patientonsetage:[65+TO+*]"  # 65 and older
}
```

#### 现场存在
检查字段是否存在：
```python
# Records that have a patient age
params = {
    "search": "_exists_:patient.patientonsetage"
}

# Records missing patient age
params = {
    "search": "_missing_:patient.patientonsetage"
}
```

### 限制参数

控制返回的结果数量（1-1000，默认 1）：
```python
params = {
    "search": "...",
    "limit": 100
}
```

**最大值**：每个请求 1000 个结果

### 跳过参数

对于分页，跳过前 N 个结果：
```python
# Get results 101-200
params = {
    "search": "...",
    "limit": 100,
    "skip": 100
}
```

**分页示例**：
```python
def get_all_results(url, search_query, api_key, max_results=5000):
    """Retrieve results with pagination."""
    all_results = []
    skip = 0
    limit = 100

    while len(all_results) < max_results:
        params = {
            "api_key": api_key,
            "search": search_query,
            "limit": limit,
            "skip": skip
        }

        response = requests.get(url, params=params)
        data = response.json()

        if "results" not in data or len(data["results"]) == 0:
            break

        all_results.extend(data["results"])

        if len(data["results"]) < limit:
            break  # No more results

        skip += limit
        time.sleep(0.25)  # Rate limiting courtesy

    return all_results[:max_results]
```

### 计数参数

按字段聚合和计算结果（而不是返回单个记录）：
```python
# Count events by country
params = {
    "search": "patient.drug.medicinalproduct:aspirin",
    "count": "occurcountry"
}
```

**回复格式**：
```json
{
  "results": [
    {"term": "us", "count": 12543},
    {"term": "ca", "count": 3421},
    {"term": "gb", "count": 2156}
  ]
}
```

#### 精确计数

添加 `.exact` 后缀以进行精确的短语计数（对于多单词字段尤其重要）：
```python
# Count exact reaction terms (not individual words)
params = {
    "search": "patient.drug.medicinalproduct:aspirin",
    "count": "patient.reaction.reactionmeddrapt.exact"
}
```

**没有 `.exact`**：计算单个单词
**使用 `.exact`**：计算完整的短语

### 排序参数

按字段对结果排序：
```python
# Sort by date, newest first
params = {
    "search": "...",
    "sort": "receivedate:desc"
}

# Sort by date, oldest first
params = {
    "search": "...",
    "sort": "receivedate:asc"
}
```

## 响应格式

### 标准响应结构

```json
{
  "meta": {
    "disclaimer": "...",
    "terms": "...",
    "license": "...",
    "last_updated": "2024-01-15",
    "results": {
      "skip": 0,
      "limit": 10,
      "total": 15234
    }
  },
  "results": [
    {
      // Individual result record
    },
    {
      // Another result record
    }
  ]
}
```

### 响应字段

- **元**：有关查询和结果的元数据
  - `disclaimer`：重要的法律免责声明
  - `terms`：使用条款 URL
  - `license`：数据许可证信息
  - `last_updated`：上次更新数据的时间
  - `results.skip`：跳过的结果数
  - `results.limit`：每页最大结果数
  - `results.total`：总匹配结果（对于大型结果集可能是近似值）

- **结果**：匹配记录数组

### 空结果

当没有结果匹配时：
```json
{
  "meta": {...},
  "results": []
}
```

### 错误响应

发生错误时：
```json
{
  "error": {
    "code": "INVALID_QUERY",
    "message": "Detailed error message"
  }
}
```

**常见错误代码**：
- `NOT_FOUND`：未找到结果 (404)
- `INVALID_QUERY`：格式错误的搜索查询 (400)
- `RATE_LIMIT_EXCEEDED`：请求太多 (429)
- `UNAUTHORIZED`：无效的 API 密钥 (401)
- `SERVER_ERROR`：内部服务器错误 (500)

## 先进技术

### 嵌套字段查询

查询嵌套对象：
```python
# Drug adverse events where serious outcome is death
params = {
    "search": "serious:1+AND+seriousnessdeath:1"
}
```
### 多字段搜索

跨多个字段搜索：
```python
# Search drug name in multiple fields
params = {
    "search": "(patient.drug.medicinalproduct:aspirin+OR+patient.drug.openfda.brand_name:aspirin)"
}
```

### 复杂布尔逻辑

组合多个运算符：
```python
# (Aspirin OR Ibuprofen) AND (Heart Attack) AND NOT (US)
params = {
    "search": "(patient.drug.medicinalproduct:aspirin+OR+patient.drug.medicinalproduct:ibuprofen)+AND+patient.reaction.reactionmeddrapt:*heart*attack*+AND+NOT+occurcountry:us"
}
```

### 使用过滤器计数

在特定子集中计数：
```python
# Count reactions for serious events only
params = {
    "search": "serious:1",
    "count": "patient.reaction.reactionmeddrapt.exact"
}
```

## 最佳实践

### 1.查询效率

**做**：
- 使用特定字段搜索
- 计数前过滤
- 尽可能使用精确匹配
- 为大型数据集实现分页

**不要**：
- 使用过于宽泛的通配符（例如，`search=*`）
- 请求比需要更多的数据
- 跳过错误处理
- 忽略速率限制

### 2. 错误处理

始终处理常见错误：
```python
def safe_api_call(url, params):
    """Safely call FDA API with comprehensive error handling."""
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.HTTPError as e:
        if response.status_code == 404:
            return {"error": "No results found"}
        elif response.status_code == 429:
            return {"error": "Rate limit exceeded"}
        elif response.status_code == 400:
            return {"error": "Invalid query"}
        else:
            return {"error": f"HTTP error: {e}"}
    except requests.exceptions.ConnectionError:
        return {"error": "Connection failed"}
    except requests.exceptions.Timeout:
        return {"error": "Request timeout"}
    except requests.exceptions.RequestException as e:
        return {"error": f"Request error: {e}"}
```

### 3.数据验证

验证和清理数据：
```python
def clean_search_term(term):
    """Clean and prepare search term."""
    # Remove special characters that break queries
    term = term.replace('"', '\\"')  # Escape quotes
    term = term.strip()
    return term

def validate_date(date_str):
    """Validate date format (YYYYMMDD)."""
    import re
    if not re.match(r'^\d{8}$', date_str):
        raise ValueError("Date must be in YYYYMMDD format")
    return date_str
```

### 4. 缓存

对经常访问的数据实施缓存：
```python
import json
from pathlib import Path
import hashlib
import time

class FDACache:
    """Simple file-based cache for FDA API responses."""

    def __init__(self, cache_dir="fda_cache", ttl=3600):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.ttl = ttl  # Time to live in seconds

    def _get_cache_key(self, url, params):
        """Generate cache key from URL and params."""
        cache_string = f"{url}_{json.dumps(params, sort_keys=True)}"
        return hashlib.md5(cache_string.encode()).hexdigest()

    def get(self, url, params):
        """Get cached response if available and not expired."""
        key = self._get_cache_key(url, params)
        cache_file = self.cache_dir / f"{key}.json"

        if cache_file.exists():
            # Check if expired
            age = time.time() - cache_file.stat().st_mtime
            if age < self.ttl:
                with open(cache_file, 'r') as f:
                    return json.load(f)

        return None

    def set(self, url, params, data):
        """Cache response data."""
        key = self._get_cache_key(url, params)
        cache_file = self.cache_dir / f"{key}.json"

        with open(cache_file, 'w') as f:
            json.dump(data, f)

# Usage
cache = FDACache(ttl=3600)  # 1 hour cache

def cached_api_call(url, params):
    """API call with caching."""
    # Check cache
    cached = cache.get(url, params)
    if cached:
        return cached

    # Make request
    response = requests.get(url, params=params)
    data = response.json()

    # Cache result
    cache.set(url, params, data)

    return data
```

### 5. 速率限制管理

跟踪和遵守速率限制：
```python
import time
from collections import deque

class RateLimiter:
    """Track and enforce rate limits."""

    def __init__(self, max_per_minute=240):
        self.max_per_minute = max_per_minute
        self.requests = deque()

    def wait_if_needed(self):
        """Wait if necessary to stay under rate limit."""
        now = time.time()

        # Remove requests older than 1 minute
        while self.requests and now - self.requests[0] > 60:
            self.requests.popleft()

        # Check if at limit
        if len(self.requests) >= self.max_per_minute:
            sleep_time = 60 - (now - self.requests[0])
            if sleep_time > 0:
                time.sleep(sleep_time)
            self.requests.popleft()

        self.requests.append(time.time())

# Usage
rate_limiter = RateLimiter(max_per_minute=240)

def rate_limited_request(url, params):
    """Make request with rate limiting."""
    rate_limiter.wait_if_needed()
    return requests.get(url, params=params)
```

## 常见查询模式

### 模式 1：基于时间的分析
```python
# Get events from last 30 days
from datetime import datetime, timedelta

end_date = datetime.now()
start_date = end_date - timedelta(days=30)

params = {
    "search": f"receivedate:[{start_date.strftime('%Y%m%d')}+TO+{end_date.strftime('%Y%m%d')}]",
    "limit": 1000
}
```

### 模式 2：Top N 分析
```python
# Get top 10 most common reactions for a drug
params = {
    "search": "patient.drug.medicinalproduct:aspirin",
    "count": "patient.reaction.reactionmeddrapt.exact",
    "limit": 10
}
```

### 模式3：比较分析
```python
# Compare two drugs
drugs = ["aspirin", "ibuprofen"]
results = {}

for drug in drugs:
    params = {
        "search": f"patient.drug.medicinalproduct:{drug}",
        "count": "patient.reaction.reactionmeddrapt.exact",
        "limit": 10
    }
    results[drug] = requests.get(url, params=params).json()
```

## 其他资源

- **openFDA 主页**：https://open.fda.gov/
- **API 文档**：https://open.fda.gov/apis/
- **交互式 API 资源管理器**：https://open.fda.gov/apis/try-the-api/
- **服务条款**：https://open.fda.gov/terms/
- **GitHub**：https://github.com/FDA/openfda
- **状态页面**：检查 API 中断和维护

## 支持

对于疑问或问题：
- **GitHub 问题**：https://github.com/FDA/openfda/issues
- **电子邮件**：open-fda@fda.hhs.gov
- **讨论论坛**：查看 GitHub 讨论