<!-- 此文件由机器翻译自 foods.md -->

# FDA 食品数据库

该参考文献涵盖了可通过 openFDA 访问的 FDA 食品相关 API 端点。

## 概述

FDA 食品数据库提供有关食品的信息，包括不良事件和执法行动。这些数据库有助于跟踪食品安全问题、召回和消费者投诉。

## 可用端点

### 1. 食品不良事件

**端点**：`https://api.fda.gov/food/event.json`

**目的**：获取食品、膳食补充剂和化妆品的不良事件报告。

**数据来源**：CAERS（CFSAN 不良事件报告系统）

**关键字段**：
- `date_started` - 不良事件开始时
- `date_created` - 创建报告时
- `report_number` - 唯一报告标识符
- `outcomes` - 事件结果（例如住院、死亡）
- `reactions` - 报告不良反应/症状
- `consumer.age` - 消费者年龄
- `consumer.age_unit` - 年龄单位（年、月等）
- `consumer.gender` - 消费者性别
- `products` - 涉及的产品数组
- `products.name_brand` - 产品品牌名称
- `products.industry_code` - 产品类别代码
- `products.industry_name` - 产品类别名称
- `products.role` - 产品角色（可疑、伴随）

**产品类别（行业名称）**：
- 烘焙产品/面团/混合物/糖霜
- 饮料（咖啡、茶、软饮料等）
- 膳食补充剂
- 冰淇淋产品
- 化妆品
- 维生素和营养补充剂
- 许多其他人

**常见用例**：
- 食品安全监测
- 膳食补充剂监测
- 不良事件趋势分析
- 产品安全评估
- 消费者投诉追踪

**查询示例**：
```python
import requests

api_key = "YOUR_API_KEY"
url = "https://api.fda.gov/food/event.json"

# Find adverse events for dietary supplements
params = {
    "api_key": api_key,
    "search": "products.industry_name:Dietary+Supplements",
    "limit": 10
}

response = requests.get(url, params=params)
data = response.json()
```

<<<代码块_1>>>

<<<代码块_2>>>

<<<代码块_3>>>

### 2. 食品执法报告

**端点**：`https://api.fda.gov/food/enforcement.json`

**目的**：获取 FDA 发布的食品召回执行报告。

**数据来源**：FDA 执法报告

**关键字段**：
- `status` - 当前状态（正在进行、已完成、已终止）
- `recall_number` - 唯一召回标识符
- `classification` - I、II 或 III 类
- `product_description` - 召回食品的描述
- `reason_for_recall` - 产品为何被召回
- `product_quantity` - 召回的产品数量
- `code_info` - 批号、批次代码、UPC
- `distribution_pattern` - 地理分布
- `recalling_firm` - 公司进行召回
- `recall_initiation_date` - 调用开始时
- `report_date` - 当 FDA 收到通知时
- `voluntary_mandated` - 自愿或 FDA 强制召回
- `city` - 回忆公司城市
- `state` - 调用公司状态
- `country` - 调用公司国家/地区
- `initial_firm_notification` - 公司如何收到通知

**分类级别**：
- **I 类**：可能导致严重健康问题或死亡的危险或有缺陷的产品（例如，具有严重风险的未申报过敏原、肉毒杆菌污染）
- **II 类**：可能导致暂时健康问题或构成轻微威胁的产品（例如轻微过敏原问题、质量缺陷）
- **III 类**：不太可能引起不良健康反应但违反 FDA 法规的产品（例如标签错误、质量问题）

**常见召回原因**：
- 未申报的过敏原（牛奶、鸡蛋、花生、坚果、大豆、小麦、鱼、贝类、芝麻）
- 微生物污染（李斯特菌、沙门氏菌、大肠杆菌等）
- 异物污染（金属、塑料、玻璃）
- 标签错误
- 加工/包装不当
- 化学污染

**常见用例**：
- 食品安全监测
- 供应链风险管理
- 过敏原追踪
- 零售商召回协调
- 消费者安全警报

**查询示例**：
<<<代码块_4>>>

<<<代码块_5>>>

<<<代码块_6>>>

```python
# Get recalls by specific company
params = {
    "api_key": api_key,
    "search": "recalling_firm:*General+Mills*",
    "limit": 20
}
```

```python
# Find ongoing recalls
params = {
    "api_key": api_key,
    "search": "status:Ongoing",
    "limit": 100
}
```

```python
# Search by product type
params = {
    "api_key": api_key,
    "search": "product_description:*ice+cream*",
    "limit": 25
}
```

## 集成技巧

### 过敏原监测系统

```python
def monitor_allergen_recalls(allergens, api_key, days_back=30):
    """
    Monitor food recalls for specific allergens.

    Args:
        allergens: List of allergens to monitor (e.g., ["peanut", "milk", "soy"])
        api_key: FDA API key
        days_back: Number of days to look back

    Returns:
        List of matching recalls
    """
    import requests
    from datetime import datetime, timedelta

    # Calculate date range
    end_date = datetime.now()
    start_date = end_date - timedelta(days=days_back)
    date_range = f"[{start_date.strftime('%Y%m%d')}+TO+{end_date.strftime('%Y%m%d')}]"

    url = "https://api.fda.gov/food/enforcement.json"
    all_recalls = []

    for allergen in allergens:
        params = {
            "api_key": api_key,
            "search": f"reason_for_recall:*{allergen}*+AND+report_date:{date_range}",
            "limit": 100
        }

        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            if "results" in data:
                for result in data["results"]:
                    result["detected_allergen"] = allergen
                    all_recalls.append(result)

    return all_recalls
```

### 不良事件分析

```python
def analyze_product_adverse_events(product_name, api_key):
    """
    Analyze adverse events for a specific food product.

    Args:
        product_name: Product name or partial name
        api_key: FDA API key

    Returns:
        Dictionary with analysis results
    """
    import requests
    from collections import Counter

    url = "https://api.fda.gov/food/event.json"
    params = {
        "api_key": api_key,
        "search": f"products.name_brand:*{product_name}*",
        "limit": 1000
    }

    response = requests.get(url, params=params)
    data = response.json()

    if "results" not in data:
        return {"error": "No results found"}

    results = data["results"]

    # Extract all reactions
    all_reactions = []
    all_outcomes = []

    for event in results:
        if "reactions" in event:
            all_reactions.extend(event["reactions"])
        if "outcomes" in event:
            all_outcomes.extend(event["outcomes"])

    # Count frequencies
    reaction_counts = Counter(all_reactions)
    outcome_counts = Counter(all_outcomes)

    return {
        "total_events": len(results),
        "most_common_reactions": reaction_counts.most_common(10),
        "outcome_distribution": dict(outcome_counts),
        "serious_outcomes": sum(1 for o in all_outcomes if o in ["Hospitalization", "Death", "Disability"])
    }
```

### 召回警报系统

```python
def get_recent_recalls_by_state(state_code, api_key, days=7):
    """
    Get recent food recalls for products distributed in a specific state.

    Args:
        state_code: Two-letter state code (e.g., "CA", "NY")
        api_key: FDA API key
        days: Number of days to look back

    Returns:
        List of recent recalls affecting the state
    """
    import requests
    from datetime import datetime, timedelta

    url = "https://api.fda.gov/food/enforcement.json"

    # Calculate date range
    end_date = datetime.now()
    start_date = end_date - timedelta(days=days)
    date_range = f"[{start_date.strftime('%Y%m%d')}+TO+{end_date.strftime('%Y%m%d')}]"

    params = {
        "api_key": api_key,
        "search": f"distribution_pattern:*{state_code}*+AND+report_date:{date_range}",
        "limit": 100,
        "sort": "report_date:desc"
    }

    response = requests.get(url, params=params)
    if response.status_code == 200:
        data = response.json()
        return data.get("results", [])
    return []
```

## 最佳实践

1. **监控过敏原召回** - 对于食品服务和零售至关重要
2. **检查分销模式** - 召回可能是地区性或全国性的
3. **跟踪召回状态** - 状态从“正在进行”更改为“已完成”
4. **按分类过滤** - 优先考虑 I 类召回以立即采取行动
5. **使用日期范围** - 关注最近发生的事件以获取操作相关性
6. **交叉引用产品** - 同一产品可能同时出现在不良事件和执法中
7. **仔细解析 code_info** - 批号和 UPC 的格式有所不同
8. **考虑产品类别** - 行业代码有助于对产品进行分类
9. **跟踪严重后果** - 住院和死亡需要立即关注
10. **实施警报系统** - 自动监控关键产品/过敏原

## 需要监测的常见过敏原

FDA 认可必须申报的 9 种主要食品过敏原：
1. 牛奶
2. 鸡蛋
3. 鱼
4.甲壳类贝类
5.坚果
6. 花生
7. 小麦
8. 大豆
9. 芝麻

这些因素占食物过敏的 90% 以上，也是 I 级召回的最常见原因。

## 其他资源

- OpenFDA 食品 API 文档：https://open.fda.gov/apis/food/
- CFSAN 不良事件报告：https://www.fda.gov/food/compliance-enforcement-food/cfsan-adverse-event-reporting-system-caers
- 食品召回：https://www.fda.gov/safety/recalls-market-withdrawals-safety-alerts
- API 基础知识：请参阅此参考目录中的 `api_basics.md`
- Python 示例：参见 `scripts/fda_food_query.py`