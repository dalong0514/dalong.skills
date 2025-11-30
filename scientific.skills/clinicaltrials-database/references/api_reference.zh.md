<!-- 此文件由机器翻译自 api_reference.md -->

# ClinicalTrials.gov API v2 参考文档

## 概述

ClinicalTrials.gov API v2 是一种现代 REST API，提供对 ClinicalTrials.gov 数据库的编程访问，该数据库包含有关在世界各地进行的临床研究的信息。该 API 遵循 OpenAPI 规范 3.0，提供 JSON 和 CSV 响应格式。

**基本网址：** `https://clinicaltrials.gov/api/v2`

**API版本：** 2.0（2024年3月发布，取代经典API）

## 身份验证和速率限制

- **身份验证：** 不需要（公共 API）
- **速率限制：** 每个 IP 地址每分钟大约 50 个请求
- **响应格式：** JSON（默认）或 CSV
- **标准：** 对日期使用 ISO 8601，对富文本使用 CommonMark Markdown

## 核心端点

### 1. 搜索研究

**端点：** `GET /api/v2/studies`

使用各种查询参数和过滤器搜索临床试验。

**查询参数：**

|参数|类型 |描述 |示例|
|------------|------|-------------|---------|
| `query.cond` |字符串|疾病或状况搜索 | `lung cancer`、`diabetes` |
| `query.intr` |字符串|治疗或干预搜索| `Pembrolizumab`、`exercise` |
| `query.locn` |字符串|地理位置过滤 | `New York`、`California, USA` |
| `query.spons` |字符串|赞助商或合作者姓名 | `National Cancer Institute` |
| `query.term` |字符串|通用全文检索| `breast cancer treatment` |
| `filter.overallStatus` |字符串|基于状态的过滤（逗号分隔）| `RECRUITING,NOT_YET_RECRUITING` |
| `filter.ids` |字符串| NCT ID 交叉过滤（逗号分隔） | `NCT04852770,NCT01728545` |
| `filter.phase` |字符串|研究阶段过滤| `PHASE1,PHASE2` |
| `sort` |字符串|结果排序 | `LastUpdatePostDate:desc` |
| `pageSize` |整数 |每页结果（最多 1000 个）| `100` |
| `pageToken` |字符串|先前响应中的分页标记 | `<token>` |
| `format` |字符串|响应格式（`json` 或 `csv`）| `json` |

**有效状态值：**
- `RECRUITING` - 目前正在招募参与者
- `NOT_YET_RECRUITING` - 尚未开放招募
- `ENROLLING_BY_INVITATION` - 仅通过邀请注册
- `ACTIVE_NOT_RECRUITING` - 活跃但不再招募
- `SUSPENDED` - 暂时停止
- `TERMINATED` - 过早停止
- `COMPLETED` - 研究已结束
- `WITHDRAWN` - 注册前撤回

**有效相位值：**
- `EARLY_PHASE1` - 早期阶段 1（以前的阶段 0）
- `PHASE1` - 第 1 阶段
- `PHASE2` - 第 2 阶段
- `PHASE3` - 第 3 阶段
- `PHASE4` - 第 4 阶段
- `NA` - 不适用

**排序选项：**
- `LastUpdatePostDate:asc` / `LastUpdatePostDate:desc` - 按上次更新日期排序
- `EnrollmentCount:asc` / `EnrollmentCount:desc` - 按注册数量排序
- `StartDate:asc` / `StartDate:desc` - 按开始日期排序
- `StudyFirstPostDate:asc` / `StudyFirstPostDate:desc` - 按首次发布日期排序

**请求示例：**
```bash
curl "https://clinicaltrials.gov/api/v2/studies?query.cond=lung+cancer&filter.overallStatus=RECRUITING&pageSize=10&format=json"
```

**响应结构示例：**
<<<代码块_1>>>

### 2. 获取研究详细信息

**端点：** `GET /api/v2/studies/{NCT_ID}`

检索有关特定临床试验的综合信息。

**路径参数：**

|参数|类型 |描述 |示例|
|------------|------|-------------|---------|
| `NCT_ID` |字符串|唯一的NCT标识符| `NCT04852770` |

**查询参数：**

|参数|类型 |描述 |示例|
|------------|------|-------------|---------|
| `format` |字符串|响应格式（`json` 或 `csv`）| `json` |

**请求示例：**
<<<代码块_2>>>

## 响应数据结构

API 返回组织成分层模块的研究数据。关键部分包括：

### 协议部分

核心研究信息和设计：

- **identificationModule** - NCT ID、正式头衔、简短头衔、组织
- **statusModule** - 总体状态、开始日期、完成日期、上次更新
- **sponsorCollaboratorsModule** - 主要发起人、合作者、责任方
- **descriptionModule** - 简要总结，详细描述
- **conditionsModule** - 正在研究的条件
- **designModule** - 研究类型、阶段、注册信息、设计详细信息
- **armsInterventionsModule** - 研究武器和干预措施
- **结果模块** - 主要和次要结果
- **资格模块** - 纳入/排除标准、年龄/性别要求
- **contactsLocationsModule** - 总体联系人、学习地点
- **referencesModule** - 参考文献、链接、引文

### 派生部分

计算/得出的信息：

- **miscInfoModule** - 版本持有者，删除的国家/地区
- **conditionBrowseModule** - 条件网格术语
- **interventionBrowseModule** - 干预网格术语

### 结果部分

研究结果（如有）：

- **participantFlowModule** - 研究中的参与者流程
- **baselineCharacteristicsModule** - 基线参与者特征
- **outcomeMeasuresModule** - 结果测量结果
- **adverseEventsModule** - 不良事件数据

### 有结果

指示结果是否可用于研究的布尔值。

## 常见用例

### 用例 1：查找针对某个条件的招聘试验

搜索目前正在招募特定疾病或病症参与者的试验：

<<<代码块_3>>>

### 用例 2：按干预/药物搜索

查找测试特定干预措施或药物的试验：

<<<代码块_4>>>

### 用例 3：地理搜索

查找特定位置的试验：

<<<代码块_5>>>

### 用例 4：检索完整的研究详细信息

获取有关特定试验的全面信息：

<<<代码块_6>>>

### 用例 5：结果分页

使用分页处理大型结果集：

```python
all_studies = []
page_token = None

while True:
    params = {
        "query.cond": "cancer",
        "pageSize": 1000
    }
    if page_token:
        params['pageToken'] = page_token

    response = requests.get("https://clinicaltrials.gov/api/v2/studies", params=params)
    data = response.json()

    all_studies.extend(data['studies'])

    # Check if there are more pages
    page_token = data.get('pageToken')
    if not page_token:
        break

print(f"Retrieved {len(all_studies)} total studies")
```

### 用例 6：导出到 CSV

检索 CSV 格式的数据进行分析：

```python
params = {
    "query.cond": "alzheimer",
    "format": "csv",
    "pageSize": 100
}

response = requests.get("https://clinicaltrials.gov/api/v2/studies", params=params)
csv_data = response.text

# Save to file
with open("alzheimer_trials.csv", "w") as f:
    f.write(csv_data)
```

## 错误处理

### 常见 HTTP 状态代码

- **200 OK** - 请求成功
- **400 错误请求** - 无效参数或格式错误的请求
- **404 未找到** - NCT ID 未找到
- **429 请求过多** - 超出速率限制
- **500 内部服务器错误** - 服务器错误

### 错误响应示例

```json
{
  "error": {
    "code": 400,
    "message": "Invalid parameter: filter.overallStatus must be one of: RECRUITING, NOT_YET_RECRUITING, ..."
  }
}
```

### 错误处理的最佳实践

```python
import requests
import time

def search_with_retry(params, max_retries=3):
    for attempt in range(max_retries):
        try:
            response = requests.get(
                "https://clinicaltrials.gov/api/v2/studies",
                params=params,
                timeout=30
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 429:
                # Rate limited - wait and retry
                wait_time = 60  # Wait 1 minute
                print(f"Rate limited. Waiting {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                raise
        except requests.exceptions.RequestException as e:
            if attempt == max_retries - 1:
                raise
            time.sleep(2 ** attempt)  # Exponential backoff

    raise Exception("Max retries exceeded")
```

## 数据标准

### 日期格式

所有日期均使用 ISO 8601 格式和结构化对象：

```json
"lastUpdatePostDateStruct": {
  "date": "2024-03-15",
  "type": "ACTUAL"
}
```

### 富文本

描述性文本字段使用 CommonMark Markdown 格式，允许结构化格式：

```json
"briefSummary": "This is a **Phase 2** study evaluating:\n\n- Safety\n- Efficacy\n- Tolerability"
```

### 枚举值

许多字段使用标准化的枚举值（例如，研究状态、阶段）而不是自由格式的文本，从而提高了数据一致性和查询可靠性。

## 从经典 API 迁移

API v2 取代了经典 API（2024 年 6 月停用）。主要改进：

1. **结构化数据** - 枚举值而不是自由文本
2. **现代标准** - ISO 8601 日期、CommonMark 降价
3. **更好的性能** - 优化的查询和分页
4. **OpenAPI Spec** - 标准 API 规范格式
5. **一致的字段** - 正确输入的数字字段

有关详细迁移指南，请参阅：https://clinicaltrials.gov/data-api/about-api/api-migration