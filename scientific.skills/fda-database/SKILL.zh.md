<!-- 此文件由机器翻译自 SKILL.md -->

---
名称： FDA 数据库
描述：“查询 openFDA API 的药物、设备、不良事件、召回、监管提交（510k、PMA）、物质识别（UNII），用于 FDA 监管数据分析和安全研究。”
---

# FDA 数据库访问

## 概述

通过 openFDA 访问全面的 FDA 监管数据，openFDA 是 FDA 为公共数据集提供开放 API 的举措。使用具有标准化接口的Python查询有关药品、医疗器械、食品、动物/兽医产品和物质的信息。

**关键能力：**
- 查询药品、器械、食品、兽药产品的不良事件
- 访问产品标签、批准和监管提交
- 监控召回和执法行动
- 查找国家药品代码 (NDC) 和物质标识符 (UNII)
- 分析设备分类和间隙（510k、PMA）
- 跟踪药品短缺和供应问题
- 研究化学结构和物质关系

## 何时使用此技能

在处理以下情况时应使用此技能：
- **药物研究**：安全概况、不良事件、标签、批准、短缺
- **医疗器械监控**：不良事件、召回、510(k) 许可、PMA 批准
- **食品安全**：召回、过敏原追踪、不良事件、膳食补充剂
- **兽医**：按物种和品种划分的动物药物不良事件
- **化学/物质数据**：UNII 查找、CAS 编号映射、分子结构
- **监管分析**：审批途径、执法行动、合规跟踪
- **药物警戒**：上市后监测、安全信号检测
- **科学研究**：药物相互作用、比较安全性、流行病学研究

## 快速入门

### 1. 基本设置

```python
from scripts.fda_query import FDAQuery

# Initialize (API key optional but recommended)
fda = FDAQuery(api_key="YOUR_API_KEY")

# Query drug adverse events
events = fda.query_drug_events("aspirin", limit=100)

# Get drug labeling
label = fda.query_drug_label("Lipitor", brand=True)

# Search device recalls
recalls = fda.query("device", "enforcement",
                   search="classification:Class+I",
                   limit=50)
```

### 2. API 密钥设置

虽然 API 无需密钥即可工作，但注册提供了更高的速率限制：
- **无密钥**：240 个请求/分钟，1,000 个/天
- **带密钥**：240 个请求/分钟，120,000 个/天

注册地址：https://open.fda.gov/apis/authentication/

设置为环境变量：
<<<代码块_1>>>

### 3. 运行示例

<<<代码块_2>>>

## FDA 数据库类别

### 毒品

访问 6 个药物相关终点，涵盖从批准到上市后监测的整个药物生命周期。

**端点：**
1. **不良事件** - 副作用、错误和治疗失败的报告
2. **产品标签** - 处方信息、警告、指示
3. **NDC目录** - 国家药品代码产品信息
4. **执行报告** - 药品召回和安全行动
5. **Drugs@FDA** - 自 1939 年以来的历史批准数据
6. **药品短缺** - 当前和已解决的供应问题

**常见用例：**
<<<代码块_3>>>

**参考：** 详细文档请参见`references/drugs.md`

### 设备

访问 9 个与设备相关的端点，涵盖医疗设备安全、审批和注册。

**端点：**
1. **不良事件** - 设备故障、受伤、死亡
2. **510(k) 许可** - 上市前通知
3. **分类** - 设备类别和风险类别
4. **执行报告** - 设备召回
5. **召回** - 详细召回信息
6. **PMA** - III 类设备的上市前批准数据
7. **注册和列表** - 制造设施数据
8. **UDI** - 唯一设备标识数据库
9. **COVID-19 血清学** - 抗体测试性能数据

**常见用例：**
<<<代码块_4>>>

**参考：** 详细文档请参见`references/devices.md`

### 食物

访问 2 个食品相关端点以进行安全监控和召回。

**端点：**
1. **不良事件** - 食品、膳食补充剂和化妆品事件
2. **执法报告** - 食品召回

**常见用例：**
<<<代码块_5>>>

**参考：** 详细文档请参见`references/foods.md`

### 动物和兽医

获取包含特定物种信息的兽药不良事件数据。

**端点：**
1. **不良事件** - 按物种、品种和产品划分的动物药物副作用

**常见用例：**
<<<代码块_6>>>

**参考：** 详细文档请参见`references/animal_veterinary.md`

### 物质及其他

通过 UNII 代码、化学结构和关系访问分子级物质数据。

**端点：**
1. **物质数据** - UNII、CAS、化学结构、关系
2. **NSDE** - 历史物质数据（遗留）

**常见用例：**
```python
# UNII to CAS mapping
substance = fda.query_substance_by_unii("R16CO5Y76E")

# Search by name
results = fda.query_substance_by_name("acetaminophen")

# Get chemical structure
structure = fda.query("other", "substance",
    search="names.name:ibuprofen+AND+substanceClass:chemical")
```

**参考：** 详细文档请参见`references/other.md`

## 常见查询模式
### 模式 1：安全概况分析

结合多个数据源创建全面的安全配置文件：

```python
def drug_safety_profile(fda, drug_name):
    """Generate complete safety profile."""

    # 1. Total adverse events
    events = fda.query_drug_events(drug_name, limit=1)
    total = events["meta"]["results"]["total"]

    # 2. Most common reactions
    reactions = fda.count_by_field(
        "drug", "event",
        search=f"patient.drug.medicinalproduct:*{drug_name}*",
        field="patient.reaction.reactionmeddrapt",
        exact=True
    )

    # 3. Serious events
    serious = fda.query("drug", "event",
        search=f"patient.drug.medicinalproduct:*{drug_name}*+AND+serious:1",
        limit=1)

    # 4. Recent recalls
    recalls = fda.query_drug_recalls(drug_name=drug_name)

    return {
        "total_events": total,
        "top_reactions": reactions["results"][:10],
        "serious_events": serious["meta"]["results"]["total"],
        "recalls": recalls["results"]
    }
```

### 模式 2：时间趋势分析

使用日期范围分析一段时间内的趋势：

```python
from datetime import datetime, timedelta

def get_monthly_trends(fda, drug_name, months=12):
    """Get monthly adverse event trends."""
    trends = []

    for i in range(months):
        end = datetime.now() - timedelta(days=30*i)
        start = end - timedelta(days=30)

        date_range = f"[{start.strftime('%Y%m%d')}+TO+{end.strftime('%Y%m%d')}]"
        search = f"patient.drug.medicinalproduct:*{drug_name}*+AND+receivedate:{date_range}"

        result = fda.query("drug", "event", search=search, limit=1)
        count = result["meta"]["results"]["total"] if "meta" in result else 0

        trends.append({
            "month": start.strftime("%Y-%m"),
            "events": count
        })

    return trends
```

### 模式3：比较分析

并排比较多个产品：

```python
def compare_drugs(fda, drug_list):
    """Compare safety profiles of multiple drugs."""
    comparison = {}

    for drug in drug_list:
        # Total events
        events = fda.query_drug_events(drug, limit=1)
        total = events["meta"]["results"]["total"] if "meta" in events else 0

        # Serious events
        serious = fda.query("drug", "event",
            search=f"patient.drug.medicinalproduct:*{drug}*+AND+serious:1",
            limit=1)
        serious_count = serious["meta"]["results"]["total"] if "meta" in serious else 0

        comparison[drug] = {
            "total_events": total,
            "serious_events": serious_count,
            "serious_rate": (serious_count/total*100) if total > 0 else 0
        }

    return comparison
```

### 模式 4：跨数据库查找

跨多个端点链接数据：

```python
def comprehensive_device_lookup(fda, device_name):
    """Look up device across all relevant databases."""

    return {
        "adverse_events": fda.query_device_events(device_name, limit=10),
        "510k_clearances": fda.query_device_510k(device_name=device_name),
        "recalls": fda.query("device", "enforcement",
                           search=f"product_description:*{device_name}*"),
        "udi_info": fda.query("device", "udi",
                            search=f"brand_name:*{device_name}*")
    }
```

## 处理结果

### 响应结构

所有 API 响应都遵循以下结构：

```python
{
    "meta": {
        "disclaimer": "...",
        "results": {
            "skip": 0,
            "limit": 100,
            "total": 15234
        }
    },
    "results": [
        # Array of result objects
    ]
}
```

### 错误处理

始终处理潜在的错误：

```python
result = fda.query_drug_events("aspirin", limit=10)

if "error" in result:
    print(f"Error: {result['error']}")
elif "results" not in result or len(result["results"]) == 0:
    print("No results found")
else:
    # Process results
    for event in result["results"]:
        # Handle event data
        pass
```

### 分页

对于大型结果集，请使用分页：

```python
# Automatic pagination
all_results = fda.query_all(
    "drug", "event",
    search="patient.drug.medicinalproduct:aspirin",
    max_results=5000
)

# Manual pagination
for skip in range(0, 1000, 100):
    batch = fda.query("drug", "event",
                     search="...",
                     limit=100,
                     skip=skip)
    # Process batch
```

## 最佳实践

### 1. 使用特定搜索

**做：**
```python
# Specific field search
search="patient.drug.medicinalproduct:aspirin"
```

**不要：**
```python
# Overly broad wildcard
search="*aspirin*"
```

### 2.实施速率限制

`FDAQuery` 类自动处理速率限制，但请注意限制：
- 每分钟 240 个请求
- 每天 120,000 个请求（使用 API 密钥）

### 3. 缓存经常访问的数据

`FDAQuery` 类包含内置缓存（默认启用）：

```python
# Caching is automatic
fda = FDAQuery(api_key=api_key, use_cache=True, cache_ttl=3600)
```

### 4. 使用精确匹配进行计数

计数/聚合时，使用 `.exact` 后缀：

```python
# Count exact phrases
fda.count_by_field("drug", "event",
                  search="...",
                  field="patient.reaction.reactionmeddrapt",
                  exact=True)  # Adds .exact automatically
```

### 5. 验证输入数据

清理并验证搜索词：

```python
def clean_drug_name(name):
    """Clean drug name for query."""
    return name.strip().replace('"', '\\"')

drug_name = clean_drug_name(user_input)
```

## API 参考

有关以下内容的详细信息：
- **身份验证和速率限制** → 请参阅`references/api_basics.md`
- **药物数据库** → 请参阅`references/drugs.md`
- **设备数据库** → 请参阅`references/devices.md`
- **食品数据库** → 请参阅`references/foods.md`
- **动物/兽医数据库** → 请参阅`references/animal_veterinary.md`
- **物质数据库** → 请参阅`references/other.md`

## 脚本

### `scripts/fda_query.py`

主查询模块的 `FDAQuery` 类提供：
- 与所有 FDA 端点的统一接口
- 自动速率限制和缓存
- 错误处理和重试逻辑
- 常见的查询模式

### `scripts/fda_examples.py`

综合示例展示：
- 药物安全概况分析
- 设备监控监控
- 食品召回追踪
- 物质查找
- 药物比较分析
- 兽药分析

运行示例：
```bash
python scripts/fda_examples.py
```

## 其他资源

- **openFDA 主页**：https://open.fda.gov/
- **API 文档**：https://open.fda.gov/apis/
- **交互式 API 资源管理器**：https://open.fda.gov/apis/try-the-api/
- **GitHub 存储库**：https://github.com/FDA/openfda
- **服务条款**：https://open.fda.gov/terms/

## 支持和故障排除

### 常见问题

**问题**：超出速率限制
- **解决方案**：使用 API 密钥、实施延迟或降低请求频率

**问题**：未找到结果
- **解决方案**：尝试更广泛的搜索词、检查拼写、使用通配符

**问题**：无效的查询语法
- **解决方案**：检查 `references/api_basics.md` 中的查询语法

**问题**：结果中缺少字段
- **解决方案**：并非所有记录都包含所有字段；始终检查字段是否存在

### 获取帮助

- **GitHub 问题**：https://github.com/FDA/openfda/issues
- **电子邮件**：open-fda@fda.hhs.gov