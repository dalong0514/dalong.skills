<!-- 此文件由机器翻译自 trademark_api.md -->

# USPTO 商标 API 参考

## 概述

USPTO 为商标数据提供两个主要 API：

1. **商标状态和文件检索（TSDR）** - 检索商标案件状态和文件
2. **商标转让检索** - 检索商标转让记录

## 1. 商标状态和文档检索 (TSDR) API

### 概述

TSDR 支持以编程方式检索商标案件状态文档和信息。

**API版本：** v1.0

**基本网址：** `https://tsdrapi.uspto.gov/ts/cd/`

### 身份验证

需要在以下位置注册 API 密钥：https://account.uspto.gov/api-manager/

在请求标头中包含 API 密钥：
```
X-Api-Key: YOUR_API_KEY
```

### 端点

#### 通过序列号获取商标状态

<<<代码块_1>>>

**示例：**
<<<代码块_2>>>

#### 通过注册号获取商标状态

<<<代码块_3>>>

### 响应格式

返回包含全面商标信息的 JSON：

<<<代码块_4>>>

### 关键数据字段

- **申请信息：**
  - `ApplicationNumber` - 序列号
  - `ApplicationDate` - 提交日期
  - `ApplicationType` - 类型（TEAS Plus、TEAS Standard 等）

- **注册信息：**
  - `RegistrationNumber` - 注册号（如果已注册）
  - `RegistrationDate` - 注册日期

- **标记信息：**
  - `MarkVerbalElementText` - 标记文本
  - `MarkCurrentStatusExternalDescriptionText` - 当前状态
  - `MarkCurrentStatusDate` - 状态日期
  - `MarkDrawingCode` - 标记类型（文字、设计等）

- **分类：**
  - `GoodsAndServices` - 具有类别的商品/服务数组

- **业主信息：**
  - `Owners` - 商标所有者/申请人数组

- **起诉历史：**
  - `ProsecutionHistoryEntry` - 起诉中的事件数组

### 常见状态值

- **已注册** - 标记已注册并处于活动状态
- **待处理** - 申请正在审查中
- **放弃** - 申请/注册被放弃
- **取消** - 注册已取消
- **暂停** - 考试暂停
- **为反对而发表** - 发表，在反对期间
- **注册和更新** - 更新注册

### Python 示例

<<<代码块_5>>>

## 2. 商标转让检索API

### 概述

从 USPTO 转让数据库检索商标转让记录。显示所有权转让和担保权益。

**API版本：** v1.4

**基本网址：** `https://assignment-api.uspto.gov/trademark/`

### 身份验证

标头中需要 API 密钥：
<<<代码块_6>>>

### 搜索方法

#### 按注册号

```
GET /v1.4/assignment/application/{registration_number}
```

#### 按序列号

```
GET /v1.4/assignment/application/{serial_number}
```

#### 按受让人姓名

```
POST /v1.4/assignment/search
```

**请求正文：**
```json
{
  "criteria": {
    "assigneeName": "Company Name"
  }
}
```

### 响应格式

返回包含分配记录的 XML：

```xml
<assignments>
  <assignment>
    <reelFrame>12345/0678</reelFrame>
    <conveyanceText>ASSIGNMENT OF ASSIGNORS INTEREST</conveyanceText>
    <recordedDate>2020-01-15</recordedDate>
    <executionDate>2020-01-10</executionDate>
    <assignors>
      <assignor>
        <name>Original Owner LLC</name>
      </assignor>
    </assignors>
    <assignees>
      <assignee>
        <name>New Owner Corporation</name>
      </assignee>
    </assignees>
  </assignment>
</assignments>
```

### 关键字段

- `reelFrame` - USPTO 卷轴和框架编号
- `conveyanceText` - 交易类型
- `recordedDate` - USPTO 记录的日期
- `executionDate` - 文档执行日期
- `assignors` - 原始所有者
- `assignees` - 新所有者
- `propertyNumbers` - 受影响的序列号/注册号

### 常见的运输类型

- **转让人权益的转让** - 所有权转让
- **担保协议** - 抵押品/担保权益
- **合并** - 公司合并
- **更改姓名** - 更改姓名
- **部分权益转让** - 部分所有权转让

### Python 示例

```python
import requests
import xml.etree.ElementTree as ET

def search_trademark_assignments(registration_number, api_key):
    """Search assignments for a trademark registration."""
    url = f"https://assignment-api.uspto.gov/trademark/v1.4/assignment/application/{registration_number}"
    headers = {"X-Api-Key": api_key}

    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.text  # Returns XML
    else:
        raise Exception(f"API error: {response.status_code}")

# Usage
xml_data = search_trademark_assignments("5678901", "YOUR_API_KEY")
root = ET.fromstring(xml_data)

for assignment in root.findall('.//assignment'):
    reel_frame = assignment.find('reelFrame').text
    recorded_date = assignment.find('recordedDate').text
    conveyance = assignment.find('conveyanceText').text

    assignor = assignment.find('.//assignor/name').text
    assignee = assignment.find('.//assignee/name').text

    print(f"{recorded_date}: {assignor} -> {assignee}")
    print(f"  Type: {conveyance}")
    print(f"  Reel/Frame: {reel_frame}\n")
```

## 用例

### 1. 监控商标状态

检查待处理申请或注册的状态：

```python
def check_trademark_health(serial_number, api_key):
    """Check if trademark needs attention."""
    data = get_trademark_status(serial_number, api_key)
    tm = data['TradeMarkAppln']

    status = tm['MarkCurrentStatusExternalDescriptionText']
    alerts = []

    if 'ABANDON' in status:
        alerts.append("⚠️ ABANDONED")
    elif 'PUBLISHED' in status:
        alerts.append("📢 In opposition period")
    elif 'SUSPENDED' in status:
        alerts.append("⏸️ Examination suspended")
    elif 'REGISTERED' in status:
        alerts.append("✅ Active")

    return alerts
```

### 2. 跟踪所有权变更

监控分配记录以了解所有权变更：

```python
def get_current_owner(registration_number, api_key):
    """Find current trademark owner from assignment records."""
    xml_data = search_trademark_assignments(registration_number, api_key)
    root = ET.fromstring(xml_data)

    assignments = []
    for assignment in root.findall('.//assignment'):
        date = assignment.find('recordedDate').text
        assignee = assignment.find('.//assignee/name').text
        assignments.append((date, assignee))

    # Most recent assignment
    if assignments:
        assignments.sort(reverse=True)
        return assignments[0][1]
    return None
```

### 3. 投资组合管理

分析商标组合：

```python
def analyze_portfolio(serial_numbers, api_key):
    """Analyze status of multiple trademarks."""
    results = {
        'active': 0,
        'pending': 0,
        'abandoned': 0,
        'expired': 0
    }

    for sn in serial_numbers:
        data = get_trademark_status(sn, api_key)
        status = data['TradeMarkAppln']['MarkCurrentStatusExternalDescriptionText']

        if 'REGISTERED' in status:
            results['active'] += 1
        elif 'PENDING' in status or 'PUBLISHED' in status:
            results['pending'] += 1
        elif 'ABANDON' in status:
            results['abandoned'] += 1
        elif 'EXPIRED' in status or 'CANCELLED' in status:
            results['expired'] += 1

    return results
```

## 速率限制和最佳实践

1. **尊重速率限制** - 使用指数退避实现重试逻辑
2. **缓存响应** - 商标数据不经常更改
3. **批处理** - 随着时间的推移分散大型投资组合的请求
4. **错误处理** - 优雅地处理丢失的数据（并非所有标记都有所有字段）
5. **数据验证** - 在 API 调用之前验证序列号/注册号

## 与其他数据集成

将商标数据与其他来源相结合：

- **TSDR + 分配** - 当前状态 + 所有权历史记录
- **多个标记** - 分析一个家族中的相关标记
- **专利数据** - 交叉引用知识产权组合

## 资源

- **TSDR API**：https://developer.uspto.gov/api-catalog/tsdr-data-api
- **分配 API**：https://developer.uspto.gov/api-catalog/trademark-assignment-search-data-api
- **API 密钥注册**：https://account.uspto.gov/api-manager/
- **商标搜索**：https://tmsearch.uspto.gov/
- **Swagger 文档**：https://developer.uspto.gov/swagger/tsdr-api-v1