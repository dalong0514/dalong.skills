<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： USPTO 数据库
描述：“访问 USPTO API 进行专利/商标检索、审查历史 (PEDS)、作业、引文、审查意见通知书、TSDR、知识产权分析和现有技术检索。”
---

# 美国专利商标局数据库

## 概述

USPTO 为专利和商标数据提供专门的 API。按关键字/发明人/受让人搜索专利，通过 PEDS 检索审查历史记录，跟踪作业，分析引文和审查意见通知书，访问商标的 TSDR，进行知识产权分析和现有技术搜索。

## 何时使用此技能

该技能应该在以下情况下使用：

- **专利搜索**：按关键字、发明人、受让人、分类或日期查找专利
- **专利详细信息**：检索完整的专利数据，包括权利要求、摘要、引文
- **商标搜索**：通过序列号或注册号查找商标
- **商标状态**：检查商标状态、所有权和申请历史
- **审查历史**：从 PEDS（专利审查数据系统）访问专利审查数据
- **审查意见书**：检索审查意见书文本、引文和驳回
- **转让**：跟踪专利/商标所有权转让
- **引文**：分析专利引文（向前和向后）
- **诉讼**：访问专利诉讼记录
- **组合分析**：分析公司或发明人的专利/商标组合

## USPTO API 生态系统

USPTO 提供多种专用 API 来满足不同的数据需求：

### 核心 API

1. **PatentSearch API** - 基于现代 ElasticSearch 的专利搜索（于 2025 年 5 月取代旧版 PatentsView）
   - 按关键字、发明人、受让人、分类、日期搜索专利
   - 2025 年 6 月 30 日之前的专利数据访问权限
   - 45 个请求/分钟的速率限制
   - **基本 URL**：`https://search.patentsview.org/api/v1/`

2. **PEDS（专利审查数据系统）** - 专利审查历史
   - 1981年至今的申请状态和交易历史
   - 审查意见通知书日期和审查活动
   - 使用 `uspto-opendata-python` Python 库
   - **替换**：PAIR 批量数据 (PBD) - 停用

3. **TSDR（商标状态和文档检索）** - 商标数据
   - 商标状态、所有权、申请历史
   - 按序列号或注册号搜索
   - **基本 URL**：`https://tsdrapi.uspto.gov/ts/cd/`

### 附加 API

4. **专利转让检索** - 所有权记录和转让
5. **商标转让检索** - 商标所有权变更
6. **丰富的引文 API** - 专利引文分析
7. **审查意见文本检索** - 审查意见全文
8. **审查意见通知书** - 审查意见通知书的引用
9. **审查意见驳回** - 驳回原因和类型
10. **PTAB API** - 专利审判和上诉委员会程序
11. **专利诉讼案件** - 联邦地方法院诉讼数据
12. **癌症登月数据集** - 癌症相关专利

## 快速入门

### API 密钥注册

所有 USPTO API 都需要 API 密钥。注册于：
**https://account.uspto.gov/api-manager/**

将 API 密钥设置为环境变量：
```bash
export USPTO_API_KEY="your_api_key_here"
```

### 帮助脚本

该技能包括常见操作的Python脚本：

- **`scripts/patent_search.py`** - 用于搜索专利的 PatentSearch API 客户端
- **`scripts/peds_client.py`>** - 用于检查历史记录的 PEDS 客户端
- **`scripts/trademark_client.py`>** - 商标数据的 TSDR 客户端

## 任务 1：检索专利

### 使用 PatentSearch API

PatentSearch API 使用带有各种运算符的 JSON 查询语言来实现灵活的搜索。

#### 基本专利检索示例

**按摘要中的关键词搜索：**
<<<代码块_1>>>

**按发明人搜索：**
<<<代码块_2>>>

**按受让人/公司搜索：**
<<<代码块_3>>>

**按日期范围搜索：**
<<<代码块_4>>>

**按CPC分类搜索：**
<<<代码块_5>>>

#### 高级专利检索

将多个条件与逻辑运算符结合起来：

<<<代码块_6>>>

#### 直接 API 使用

对于复杂的查询，直接使用API：

```python
import requests

url = "https://search.patentsview.org/api/v1/patent"
headers = {
    "X-Api-Key": "YOUR_API_KEY",
    "Content-Type": "application/json"
}

query = {
    "q": {
        "_and": [
            {"patent_date": {"_gte": "2024-01-01"}},
            {"assignee_organization": {"_text_any": ["Google", "Alphabet"]}},
            {"cpc_subclass_id": ["G06N", "H04N"]}
        ]
    },
    "f": ["patent_number", "patent_title", "patent_date", "inventor_name"],
    "s": [{"patent_date": "desc"}],
    "o": {"per_page": 100, "page": 1}
}

response = requests.post(url, headers=headers, json=query)
results = response.json()
```

### 查询运算符

- **平等**：`{"field": "value"}` 或 `{"field": {"_eq": "value"}}`
- **比较**：`_gt`、`_gte`、`_lt`、`_lte`、`_neq`
- **文本搜索**：`_text_all`、`_text_any`、`_text_phrase`
- **字符串匹配**：`_begins`、`_contains`
- **逻辑**：`_and`、`_or`、`_not`
**最佳实践**：对文本字段使用 `_text_*` 运算符（比 `_contains` 或 `_begins` 性能更高）

### 可用的专利端点

- `/patent` - 已授予专利
- `/publication` - 预授权出版物
- `/inventor` - 发明人信息
- `/assignee` - 受让人信息
- `/cpc_subclass`、`/cpc_at_issue` - CPC 分类
- `/uspc` - 美国专利分类
- `/ipc` - 国际专利分类
- `/claims`、`/brief_summary_text`、`/detail_description_text` - 文本数据（测试版）

### 参考文档

请参阅 `references/patentsearch_api.md` 以获取完整的 PatentSearch API 文档，包括：
- 所有可用端点
- 完整的现场参考
- 查询语法和示例
- 响应格式
- 速率限制和最佳实践

## 任务2：检索专利审查数据

### 使用 PEDS（专利审查数据系统）

PEDS 提供全面的起诉历史记录，包括交易事件、状态变化和审查时间表。

####安装

```bash
uv pip install uspto-opendata-python
```

#### PEDS 基本用法

**获取应用程序数据：**
```python
from scripts.peds_client import PEDSHelper

helper = PEDSHelper()

# By application number
app_data = helper.get_application("16123456")
print(f"Title: {app_data['title']}")
print(f"Status: {app_data['app_status']}")

# By patent number
patent_data = helper.get_patent("11234567")
```

**获取交易历史：**
```python
transactions = helper.get_transaction_history("16123456")

for trans in transactions:
    print(f"{trans['date']}: {trans['code']} - {trans['description']}")
```

**获取审查意见通知书：**
```python
office_actions = helper.get_office_actions("16123456")

for oa in office_actions:
    if oa['code'] == 'CTNF':
        print(f"Non-final rejection: {oa['date']}")
    elif oa['code'] == 'CTFR':
        print(f"Final rejection: {oa['date']}")
    elif oa['code'] == 'NOA':
        print(f"Notice of allowance: {oa['date']}")
```

**获取状态摘要：**
```python
summary = helper.get_status_summary("16123456")

print(f"Current status: {summary['current_status']}")
print(f"Filing date: {summary['filing_date']}")
print(f"Pendency: {summary['pendency_days']} days")

if summary['is_patented']:
    print(f"Patent number: {summary['patent_number']}")
    print(f"Issue date: {summary['issue_date']}")
```

#### 检方分析

分析起诉模式：

```python
analysis = helper.analyze_prosecution("16123456")

print(f"Total office actions: {analysis['total_office_actions']}")
print(f"Non-final rejections: {analysis['non_final_rejections']}")
print(f"Final rejections: {analysis['final_rejections']}")
print(f"Allowed: {analysis['allowance']}")
print(f"Responses filed: {analysis['responses']}")
```

### 常用交易代码

- **CTNF** - 邮寄非最终拒绝邮件
- **CTFR** - 最终拒绝邮寄
- **NOA** - 邮寄津贴通知
- **书面** - 已提交回复
- **ISS.FEE** - 发行费用支付
- **ABND** - 申请被放弃
- **AOPF** - 邮寄办公室行动

### 参考文档

请参阅 `references/peds_api.md` 以获取完整的 PEDS 文档，包括：
- 所有可用的数据字段
- 交易代码参考
- Python库的使用
- 投资组合分析示例

## 任务 3：检索和监控商标

### 使用 TSDR（商标状态和文档检索）

查看商标状态、所有权和起诉历史记录。

#### 商标基本使用

**通过序列号获取商标：**
```python
from scripts.trademark_client import TrademarkClient

client = TrademarkClient()

# By serial number
tm_data = client.get_trademark_by_serial("87654321")

# By registration number
tm_data = client.get_trademark_by_registration("5678901")
```

**获取商标状态：**
```python
status = client.get_trademark_status("87654321")

print(f"Mark: {status['mark_text']}")
print(f"Status: {status['status']}")
print(f"Filing date: {status['filing_date']}")

if status['is_registered']:
    print(f"Registration #: {status['registration_number']}")
    print(f"Registration date: {status['registration_date']}")
```

**检查商标健康状况：**
```python
health = client.check_trademark_health("87654321")

print(f"Mark: {health['mark']}")
print(f"Status: {health['status']}")

for alert in health['alerts']:
    print(alert)

if health['needs_attention']:
    print("⚠️  This mark needs attention!")
```

#### 商标组合监控

监控多个商标：

```python
def monitor_portfolio(serial_numbers, api_key):
    """Monitor trademark portfolio health."""
    client = TrademarkClient(api_key)

    results = {
        'active': [],
        'pending': [],
        'problems': []
    }

    for sn in serial_numbers:
        health = client.check_trademark_health(sn)

        if 'REGISTERED' in health['status']:
            results['active'].append(health)
        elif 'PENDING' in health['status'] or 'PUBLISHED' in health['status']:
            results['pending'].append(health)
        elif health['needs_attention']:
            results['problems'].append(health)

    return results
```

### 常见商标状态

- **已注册** - 有效注册商标
- **待定** - 正在审查中
- **发表反对意见** - 在反对期间
- **放弃** - 申请被放弃
- **取消** - 注册已取消
- **暂停** - 考试暂停
- **注册和更新** - 更新注册

### 参考文档

请参阅 `references/trademark_api.md` 了解完整的商标 API 文档，包括：
- TSDR API 参考
- 商标转让搜索API
- 所有状态代码
- 起诉历史记录访问
- 所有权追踪

## 任务 4：跟踪分配和所有权

### 专利和商标转让

专利和商标都具有用于跟踪所有权变更的转让搜索 API。

#### 专利转让 API

**基本 URL**：`https://assignment-api.uspto.gov/patent/v1.4/`

**按专利号搜索：**
```python
import requests
import xml.etree.ElementTree as ET

def get_patent_assignments(patent_number, api_key):
    url = f"https://assignment-api.uspto.gov/patent/v1.4/assignment/patent/{patent_number}"
    headers = {"X-Api-Key": api_key}

    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.text  # Returns XML

assignments_xml = get_patent_assignments("11234567", api_key)
root = ET.fromstring(assignments_xml)

for assignment in root.findall('.//assignment'):
    recorded_date = assignment.find('recordedDate').text
    assignor = assignment.find('.//assignor/name').text
    assignee = assignment.find('.//assignee/name').text
    conveyance = assignment.find('conveyanceText').text

    print(f"{recorded_date}: {assignor} → {assignee}")
    print(f"  Type: {conveyance}\n")
```

**按公司名称搜索：**
```python
def find_company_patents(company_name, api_key):
    url = "https://assignment-api.uspto.gov/patent/v1.4/assignment/search"
    headers = {"X-Api-Key": api_key}
    data = {"criteria": {"assigneeName": company_name}}

    response = requests.post(url, headers=headers, json=data)
    return response.text
```

### 常见作业类型

- **转让人权益的转让** - 所有权转让
- **担保协议** - 抵押品/担保权益
- **合并** - 公司合并
- **更改姓名** - 更改姓名
- **部分权益的转让** - 部分所有权

## 任务 5：访问其他 USPTO 数据

### 审查意见通知书、传票和诉讼

多个专用 API 提供额外的专利数据。

#### 办公通知文本检索

使用申请号检索审查意见通知全文。与 PEDS 集成以识别存在哪些审查意见，然后检索全文。

#### 丰富的引文 API

分析专利引用：
- 前向引用（引用本专利的专利）
- 向后引用（引用的现有技术）
- 审查员与申请人的引文
- 引文上下文

####专利诉讼案件API

访问联邦地区法院专利诉讼记录：
- 74,623+ 诉讼记录
- 已主张的专利
- 派对和场地
- 案件结果

#### PTAB API

专利审判和上诉委员会程序：
- 多方复审（IPR）
- 授权后审查（PGR）
- 上诉决定

### 参考文档

有关以下内容的综合文档，请参阅 `references/additional_apis.md`：
- 丰富的引文API
- 审查意见通知 API（文本、引文、驳回）
- 专利诉讼案件API
- PTAB API
- 癌症登月数据集
- OCE 状态/事件代码

## 完整分析示例

### 综合专利分析

结合多个 API 以获得完整的专利情报：

```python
def comprehensive_patent_analysis(patent_number, api_key):
    """
    Full patent analysis using multiple USPTO APIs.
    """
    from scripts.patent_search import PatentSearchClient
    from scripts.peds_client import PEDSHelper

    results = {}

    # 1. Get patent details
    patent_client = PatentSearchClient(api_key)
    patent_data = patent_client.get_patent(patent_number)
    results['patent'] = patent_data

    # 2. Get examination history
    peds = PEDSHelper()
    results['prosecution'] = peds.analyze_prosecution(patent_number)
    results['status'] = peds.get_status_summary(patent_number)

    # 3. Get assignment history
    import requests
    assign_url = f"https://assignment-api.uspto.gov/patent/v1.4/assignment/patent/{patent_number}"
    assign_resp = requests.get(assign_url, headers={"X-Api-Key": api_key})
    results['assignments'] = assign_resp.text if assign_resp.status_code == 200 else None

    # 4. Analyze results
    print(f"\n=== Patent {patent_number} Analysis ===\n")
    print(f"Title: {patent_data['patent_title']}")
    print(f"Assignee: {', '.join(patent_data.get('assignee_organization', []))}")
    print(f"Issue Date: {patent_data['patent_date']}")

    print(f"\nProsecution:")
    print(f"  Office Actions: {results['prosecution']['total_office_actions']}")
    print(f"  Rejections: {results['prosecution']['non_final_rejections']} non-final, {results['prosecution']['final_rejections']} final")
    print(f"  Pendency: {results['prosecution']['pendency_days']} days")

    # Analyze citations
    if 'cited_patent_number' in patent_data:
        print(f"\nCitations:")
        print(f"  Cites: {len(patent_data['cited_patent_number'])} patents")
    if 'citedby_patent_number' in patent_data:
        print(f"  Cited by: {len(patent_data['citedby_patent_number'])} patents")

    return results
```

## 最佳实践

1. **API密钥管理**
   - 将 API 密钥存储在环境变量中
   - 切勿将密钥提交给版本控制
   - 在所有 USPTO API 中使用相同的密钥

2. **速率限制**
   - 专利搜索：45 次请求/分钟
   - 对速率限制错误实施指数退避
   - 尽可能缓存响应

3. **查询优化**
   - 对文本字段使用 `_text_*` 运算符（性能更高）
   - 仅请求所需的字段以减少响应大小
   - 使用日期范围缩小搜索范围

4. **数据处理**
   - 并非所有字段都填充了所有专利/商标
   - 优雅地处理丢失的数据
   - 一致地解析日期

5. **组合API**
   - 使用 PatentSearch 进行发现
   - 使用 PEDS 获取起诉详情
   - 使用分配 API 进行所有权跟踪
   - 结合数据进行综合分析

## 重要提示

- **旧版 API 日落**：PatentsView 旧版 API 于 2025 年 5 月 1 日停止使用 - 使用 PatentSearch API
- **PAIR 批量数据停用**：使用 PEDS 代替
- **数据覆盖范围**：PatentSearch 拥有截至 2025 年 6 月 30 日的数据； 1981 年至今的 PEDS
- **文本端点**：声明和描述端点处于测试阶段，正在进行回填
- **速率限制**：遵守速率限制以避免服务中断

## 资源

### API 文档
- **专利搜索 API**：https://search.patentsview.org/docs/
- **美国专利商标局开发者门户**：https://developer.uspto.gov/
- **美国专利商标局开放数据门户**：https://data.uspto.gov/
- **API 密钥注册**：https://account.uspto.gov/api-manager/

### Python 库
- **uspto-opendata-python**：https://pypi.org/project/uspto-opendata-python/
- **美国专利商标局文件**：https://docs.ip-tools.org/uspto-opendata-python/

### 参考文件
- `references/patentsearch_api.md` - 完整的 PatentSearch API 参考
- `references/peds_api.md` - PEDS API 和库文档
- `references/trademark_api.md` - 商标 API（TSDR 和分配）
- `references/additional_apis.md` - 传票、审查意见通知书、诉讼、PTAB

### 脚本
- `scripts/patent_search.py` - PatentSearch API 客户端
- `scripts/peds_client.py` - PEDS 检查数据客户端
- `scripts/trademark_client.py` - 商标搜索客户端