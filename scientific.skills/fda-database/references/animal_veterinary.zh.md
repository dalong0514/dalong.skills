<!-- 此文件由机器翻译自 animal_veterinary.md -->

# FDA 动物和兽医数据库

该参考文献涵盖可通过 openFDA 访问的 FDA 动物和兽医 API 端点。

## 概述

FDA 动物和兽医数据库提供有关动物药物和兽药产品相关不良事件的信息。这些数据库有助于监控伴侣动物、牲畜和其他动物所用产品的安全性。

## 可用端点

### 动物药品不良事件

**端点**：`https://api.fda.gov/animalandveterinary/event.json`

**目的**：获取与动物药物相关的副作用、产品使用错误、产品质量问题和治疗失败的报告。

**数据来源**：FDA 兽医中心 (CVM) 不良事件报告系统

**关键字段**：
- `unique_aer_id_number` - 唯一的不良事件报告标识符
- `report_id` - 报告 ID 号
- `receiver.organization` - 接收报告的组织
- `receiver.street_address` - 接收者地址
- `receiver.city` - 接收城市
- `receiver.state` - 接收器状态
- `receiver.postal_code` - 收件人邮政编码
- `receiver.country` - 接收方国家/地区
- `primary_reporter` - 主要报告者类型（例如兽医、主人）
- `onset_date` - 不良事件开始日期
- `animal.species` - 受影响的动物种类
- `animal.gender` - 动物性别
- `animal.age.min` - 最低年龄
- `animal.age.max` - 最大年龄
- `animal.age.unit` - 年龄单位（天、月、年）
- `animal.age.qualifier` - 年龄限定符
- `animal.breed.is_crossbred` - 是否杂交
- `animal.breed.breed_component` - 品种
- `animal.weight.min` - 最小重量
- `animal.weight.max` - 最大重量
- `animal.weight.unit` - 重量单位
- `animal.female_animal_physiological_status` - 繁殖状态
- `animal.reproductive_status` - 绝育/绝育状态
- `drug` - 涉及的药物数组
- `drug.active_ingredients` - 活性成分
- `drug.active_ingredients.name` - 成分名称
- `drug.active_ingredients.dose` - 剂量信息
- `drug.brand_name` - 品牌名称
- `drug.manufacturer.name` - 制造商
- `drug.administered_by` - 谁服用了药物
- `drug.route` - 给药途径
- `drug.dosage_form` - 剂型
- `drug.atc_vet_code` - ATC 兽医代码
- `reaction` - 一系列不良反应
- `reaction.veddra_version` - VeDDRA 字典版本
- `reaction.veddra_term_code` - VeDDRA 术语代码
- `reaction.veddra_term_name` - VeDDRA 术语名称
- `reaction.accuracy` - 诊断准确性
- `reaction.number_of_animals_affected` - 受影响的数量
- `reaction.number_of_animals_treated` - 处理的数字
- `outcome.medical_status` - 医疗结果
- `outcome.number_of_animals_affected` - 受结果影响的动物
- `serious_ae` - 是否存在严重不良事件
- `health_assessment_prior_to_exposure.assessed_by` - 谁评估了健康状况
- `health_assessment_prior_to_exposure.condition` - 健康状况
- `treated_for_ae` - 是否已处理
- `time_between_exposure_and_onset` - 开始时间
- `duration.unit` - 持续时间单位
- `duration.value` - 持续时间值

**常见动物种类**：
- 狗（犬狼疮）
- 猫（Felis catus）
- 马（Equus caballus）
- 牛（Bos taurus）
- 猪（Sus scrofa Domesticus）
- 鸡（鸡内金）
- 绵羊（绵羊）
- 山羊（Capra aegagrus hircus）
- 还有许多其他人

**常见用例**：
- 兽医药物警戒
- 产品安全监控
- 不良事件趋势分析
- 药品安全性比较
- 特定物种的安全性研究
- 品种倾向研究

**查询示例**：
```python
import requests

api_key = "YOUR_API_KEY"
url = "https://api.fda.gov/animalandveterinary/event.json"

# Find adverse events in dogs
params = {
    "api_key": api_key,
    "search": "animal.species:Dog",
    "limit": 10
}

response = requests.get(url, params=params)
data = response.json()
```

<<<代码块_1>>>

<<<代码块_2>>>

<<<代码块_3>>>

<<<代码块_4>>>

<<<代码块_5>>>

<<<代码块_6>>>

## VeDDRA - 药物相关事务兽医词典

药物相关事务兽医词典 (VeDDRA) 是用于不良事件报告的标准化国际兽医术语。它提供：

- 兽医不良事件的标准化术语
- 术语的层次结构
- 物种特定术语
- 国际协调

**VeDDRA 期限结构**：
- 术语按层次结构组织
- 每个术语都有一个唯一的代码
- 术语适合物种
- 存在多个版本（检查`veddra_version`字段）

## 集成技巧

### 物种特异性不良事件分析

```python
def analyze_species_adverse_events(species, drug_name, api_key):
    """
    Analyze adverse events for a specific drug in a particular species.

    Args:
        species: Animal species (e.g., "Dog", "Cat", "Horse")
        drug_name: Drug brand name or active ingredient
        api_key: FDA API key

    Returns:
        Dictionary with analysis results
    """
    import requests
    from collections import Counter

    url = "https://api.fda.gov/animalandveterinary/event.json"
    params = {
        "api_key": api_key,
        "search": f"animal.species:{species}+AND+drug.brand_name:*{drug_name}*",
        "limit": 1000
    }

    response = requests.get(url, params=params)
    data = response.json()

    if "results" not in data:
        return {"error": "No results found"}

    results = data["results"]

    # Collect reactions and outcomes
    reactions = []
    outcomes = []
    serious_count = 0

    for event in results:
        if "reaction" in event:
            for reaction in event["reaction"]:
                if "veddra_term_name" in reaction:
                    reactions.append(reaction["veddra_term_name"])

        if "outcome" in event:
            for outcome in event["outcome"]:
                if "medical_status" in outcome:
                    outcomes.append(outcome["medical_status"])

        if event.get("serious_ae") == "true":
            serious_count += 1

    reaction_counts = Counter(reactions)
    outcome_counts = Counter(outcomes)

    return {
        "total_events": len(results),
        "serious_events": serious_count,
        "most_common_reactions": reaction_counts.most_common(10),
        "outcome_distribution": dict(outcome_counts),
        "serious_percentage": round((serious_count / len(results)) * 100, 2) if len(results) > 0 else 0
    }
```

### 品种倾向研究

```python
def analyze_breed_predisposition(reaction_term, api_key, min_events=5):
    """
    Identify breed predispositions for specific adverse reactions.

    Args:
        reaction_term: VeDDRA reaction term to analyze
        api_key: FDA API key
        min_events: Minimum number of events to include breed

    Returns:
        List of breeds with event counts
    """
    import requests
    from collections import Counter

    url = "https://api.fda.gov/animalandveterinary/event.json"
    params = {
        "api_key": api_key,
        "search": f"reaction.veddra_term_name:*{reaction_term}*",
        "limit": 1000
    }

    response = requests.get(url, params=params)
    data = response.json()

    if "results" not in data:
        return []

    breeds = []
    for event in data["results"]:
        if "animal" in event and "breed" in event["animal"]:
            breed_info = event["animal"]["breed"]
            if "breed_component" in breed_info:
                if isinstance(breed_info["breed_component"], list):
                    breeds.extend(breed_info["breed_component"])
                else:
                    breeds.append(breed_info["breed_component"])

    breed_counts = Counter(breeds)

    # Filter by minimum events
    filtered_breeds = [
        {"breed": breed, "count": count}
        for breed, count in breed_counts.most_common()
        if count >= min_events
    ]

    return filtered_breeds
```

### 药物安全性比较

```python
def compare_drug_safety(drug_list, species, api_key):
    """
    Compare safety profiles of multiple drugs for a specific species.

    Args:
        drug_list: List of drug names to compare
        species: Animal species
        api_key: FDA API key

    Returns:
        Dictionary comparing drugs
    """
    import requests

    url = "https://api.fda.gov/animalandveterinary/event.json"
    comparison = {}

    for drug in drug_list:
        params = {
            "api_key": api_key,
            "search": f"animal.species:{species}+AND+drug.brand_name:*{drug}*",
            "limit": 1000
        }

        response = requests.get(url, params=params)
        data = response.json()

        if "results" in data:
            results = data["results"]
            serious = sum(1 for r in results if r.get("serious_ae") == "true")
            deaths = sum(
                1 for r in results
                if "outcome" in r
                and any(o.get("medical_status") == "Died" for o in r["outcome"])
            )

            comparison[drug] = {
                "total_events": len(results),
                "serious_events": serious,
                "deaths": deaths,
                "serious_rate": round((serious / len(results)) * 100, 2) if len(results) > 0 else 0,
                "death_rate": round((deaths / len(results)) * 100, 2) if len(results) > 0 else 0
            }

    return comparison
```

## 最佳实践

1. **使用标准物种名称** - 完整的学名或通用名效果最好
2. **考虑品种差异** - 拼写和命名可能会有所不同
3. **检查 VeDDRA 版本** - 条款可能因版本而异
4. **考虑记者偏见** - 兽医与主人的报告不同
5. **按严重事件过滤** - 关注临床显着反应
6. **考虑动物人口统计学** - 年龄、体重和繁殖状况很重要
7. **跟踪时间模式** - 可能存在季节性变化
8. **交叉参考产品** - 相同的活性成分可能有多个品牌
9. **按途径分析** - 局部给药与全身给药会影响安全性
10. **考虑物种差异** - 药物对物种的影响不同

## 报告来源

兽药不良事件报告来自：
- **兽医** - 专业医学观察
- **动物主人** - 直接观察和关注
- **制药公司** - 所需的上市后监督
- **FDA 现场工作人员** - 官方调查
- **研究机构** - 临床研究
- **其他来源** - 各不相同

不同的来源可能有不同的报告阈值和详细程度。

## 其他资源

- OpenFDA 动物和兽医 API：https://open.fda.gov/apis/animalandveterinary/
- FDA 兽医中心：https://www.fda.gov/animal-veterinary
- VeDDRA：https://www.veddra.org/
- API 基础知识：请参阅此参考目录中的 `api_basics.md`
- Python 示例：参见 `scripts/fda_animal_query.py`