<!-- 此文件由机器翻译自 devices.md -->

# FDA 医疗器械数据库

此参考涵盖通过 openFDA 访问的所有 FDA 医疗设备相关 API 端点。

## 概述

FDA 设备数据库提供有关医疗设备的信息，包括不良事件、召回、批准、注册和分类数据。医疗设备的范围从压舌板等简单物品到起搏器和手术机器人等复杂仪器。

## 设备分类系统

医疗器械根据风险分为三类：

- **I 类**：低风险（例如绷带、检查手套）
- **II 类**：中等风险（例如电动轮椅、输液泵）
- **III 类**：高风险（例如心脏瓣膜、植入式起搏器）

## 可用端点

### 1. 器械不良事件

**端点**：`https://api.fda.gov/device/event.json`

**目的**：访问记录严重伤害、死亡、故障以及医疗设备使用造成的其他不良影响的报告。

**数据来源**：制造商和用户设施设备体验 (MAUDE) 数据库

**关键字段**：
- `device.brand_name` - 设备的品牌名称
- `device.generic_name` - 通用设备名称
- `device.manufacturer_d_name` - 制造商名称
- `device.device_class` - 设备类别（1、2 或 3）
- `event_type` - 事件类型（死亡、受伤、故障、其他）
- `date_received` - FDA 收到报告的日期
- `mdr_report_key` - 唯一报告标识符
- `adverse_event_flag` - 是否报告为不良事件
- `product_problem_flag` - 是否报告产品问题
- `patient.patient_problems` - 患者问题/并发症
- `device.openfda.device_name` - 官方设备名称
- `device.openfda.medical_specialty_description` - 医学专业
- `remedial_action` - 采取的措施（召回、修理、更换等）

**常见用例**：
- 上市后监督
- 安全信号检测
- 设备比较研究
- 风险分析
- 品质提升

**查询示例**：
```python
import requests

api_key = "YOUR_API_KEY"
url = "https://api.fda.gov/device/event.json"

# Find adverse events for a specific device
params = {
    "api_key": api_key,
    "search": "device.brand_name:pacemaker",
    "limit": 10
}

response = requests.get(url, params=params)
data = response.json()
```

<<<代码块_1>>>

<<<代码块_2>>>

### 2. 设备 510(k) 间隙

**端点**：`https://api.fda.gov/device/510k.json`

**目的**：访问 510(k) 上市前通知数据，证明设备与合法销售的谓词设备等效。

**数据来源**：510(k) 上市前通知

**关键字段**：
- `k_number` - 510(k) 数字（唯一标识符）
- `applicant` - 公司提交 510(k)
- `device_name` - 设备名称
- `device_class` - 设备分类（1、2 或 3）
- `decision_date` - FDA 决定的日期
- `decision_description` - 基本等效 (SE) 或非 SE
- `product_code` - FDA 产品代码
- `statement_or_summary` - 提供的摘要类型
- `clearance_type` - 传统、特殊、缩写等。
- `expedited_review_flag` - 是否加急审核
- `advisory_committee` - 咨询委员会名称
- `openfda.device_name` - 官方设备名称
- `openfda.device_class` - 设备类描述
- `openfda.medical_specialty_description` - 医学专业
- `openfda.regulation_number` - CFR 法规编号

**常见用例**：
- 监管途径研究
- 谓词设备识别
- 市场进入分析
- 竞争情报
- 设备开发规划

**查询示例**：
<<<代码块_3>>>

<<<代码块_4>>>

<<<代码块_5>>>

### 3. 设备分类

**端点**：`https://api.fda.gov/device/classification.json`

**目的**：访问包含医疗器械名称、产品代码、医疗专业面板和分类信息的器械分类数据库。

**数据来源**：FDA 器械分类数据库

**关键字段**：
- `product_code` - 三字母 FDA 产品代码
- `device_name` - 官方设备名称
- `device_class` - 类别（1、2 或 3）
- `medical_specialty` - 医学专业（例如放射学、心血管）
- `medical_specialty_description` - 完整的专业描述
- `regulation_number` - CFR 法规编号（例如 21 CFR 870.2300）
- `review_panel` - FDA 审查小组
- `definition` - 官方设备定义
- `physical_state` - 固体、液体、气体
- `technical_method` - 操作方法
- `target_area` - 目标身体区域/系统
- `gmp_exempt_flag` - 是否免除良好生产规范
- `implant_flag` - 设备是否植入
- `life_sustain_support_flag` - 是否维持生命/支持生命

**常见用例**：
- 设备识别
- 监管要求确定
- 产品代码查找
- 分类研究
- 设备分类

**查询示例**：
<<<代码块_6>>>
```python
# Find all cardiovascular devices
params = {
    "api_key": api_key,
    "search": "medical_specialty:CV",
    "limit": 100
}
```

```python
# Get all implantable Class III devices
params = {
    "api_key": api_key,
    "search": "device_class:3+AND+implant_flag:Y",
    "limit": 50
}
```

### 4. 设备召回执行报告

**端点**：`https://api.fda.gov/device/enforcement.json`

**目的**：获取医疗器械产品召回执行报告。

**数据来源**：FDA 执法报告

**关键字段**：
- `status` - 当前状态（正在进行、已完成、已终止）
- `recall_number` - 唯一召回标识符
- `classification` - I、II 或 III 类
- `product_description` - 召回设备的描述
- `reason_for_recall` - 设备为何被召回
- `product_quantity` - 召回的产品数量
- `code_info` - 批号、序列号、型号
- `distribution_pattern` - 地理分布
- `recalling_firm` - 公司进行召回
- `recall_initiation_date` - 召回开始时
- `report_date` - 当 FDA 收到通知时
- `product_res_number` - 产品问题编号

**常见用例**：
- 质量监控
- 供应链风险管理
- 患者安全追踪
- 监管合规性
- 设备监控

**查询示例**：
```python
# Find all Class I device recalls (most serious)
params = {
    "api_key": api_key,
    "search": "classification:Class+I",
    "limit": 20,
    "sort": "report_date:desc"
}

response = requests.get("https://api.fda.gov/device/enforcement.json", params=params)
```

```python
# Search recalls by manufacturer
params = {
    "api_key": api_key,
    "search": "recalling_firm:*Philips*",
    "limit": 50
}
```

### 5. 设备召回

**端点**：`https://api.fda.gov/device/recall.json`

**目的**：获取有关设备召回的信息，以解决违反 FDA 法律或构成健康风险的问题。

**数据来源**：FDA 召回数据库

**关键字段**：
- `res_event_number` - 调用事件编号
- `product_code` - FDA 产品代码
- `openfda.device_name` - 设备名称
- `openfda.device_class` - 设备类
- `product_res_number` - 产品召回编号
- `firm_fei_number` - 公司机构标识符
- `k_numbers` - 关联的 510(k) 号码
- `pma_numbers` - 关联的 PMA 编号
- `root_cause_description` - 问题的根本原因

**常见用例**：
- 召回跟踪
- 质量调查
- 根本原因分析
- 趋势识别

**查询示例**：
```python
# Search recalls by product code
params = {
    "api_key": api_key,
    "search": "product_code:DQY",
    "limit": 20
}

response = requests.get("https://api.fda.gov/device/recall.json", params=params)
```

### 6. 上市前批准 (PMA)

**端点**：`https://api.fda.gov/device/pma.json`

**目的**：获取 FDA III 类医疗器械上市前审批流程的数据。

**数据来源**：PMA数据库

**关键字段**：
- `pma_number` - PMA 申请号（例如 P850005）
- `supplement_number` - 补充编号（如果适用）
- `applicant` - 公司名称
- `trade_name` - 商品/品牌名称
- `generic_name` - 通用名称
- `product_code` - FDA 产品代码
- `decision_date` - FDA 决定的日期
- `decision_code` - 批准状态（APPR = 已批准）
- `advisory_committee` - 咨询委员会
- `openfda.device_name` - 官方设备名称
- `openfda.device_class` - 设备类
- `openfda.medical_specialty_description` - 医学专业
- `openfda.regulation_number` - 法规编号

**常见用例**：
- 高风险设备研究
- 批准时间表分析
- 监管策略
- 市场情报
- 临床试验计划

**查询示例**：
```python
# Find PMA approvals by company
params = {
    "api_key": api_key,
    "search": "applicant:Boston+Scientific",
    "limit": 50
}

response = requests.get("https://api.fda.gov/device/pma.json", params=params)
```

```python
# Search for specific device PMAs
params = {
    "api_key": api_key,
    "search": "generic_name:*cardiac+pacemaker*",
    "limit": 10
}
```

### 7. 注册和列表

**端点**：`https://api.fda.gov/device/registrationlisting.json`

**目的**：访问医疗设备机构及其制造的设备的位置数据。

**数据源**：设备注册和列表数据库

**关键字段**：
- `registration.fei_number` - 设施建立标识符
- `registration.name` - 设施名称
- `registration.registration_number` - 注册号
- `registration.reg_expiry_date_year` - 注册到期年份
- `registration.address_line_1` - 街道地址
- `registration.city` - 城市
- `registration.state_code` - 州/省
- `registration.iso_country_code` - 国家/地区代码
- `registration.zip_code` - 邮政编码
- `products.product_code` - 设备产品代码
- `products.created_date` - 列出设备时
- `products.openfda.device_name` - 设备名称
- `products.openfda.device_class` - 设备类
- `proprietary_name` - 专有/品牌名称
- `establishment_type` - 操作类型（制造商等）

**常见用例**：
- 制造商标识
- 设施位置查找
- 供应链映射
- 尽职调查研究
- 市场分析

**查询示例**：
```python
# Find registered facilities by country
params = {
    "api_key": api_key,
    "search": "registration.iso_country_code:US",
    "limit": 100
}

response = requests.get("https://api.fda.gov/device/registrationlisting.json", params=params)
```

```python
# Search by facility name
params = {
    "api_key": api_key,
    "search": "registration.name:*Johnson*",
    "limit": 10
}
```

### 8. 唯一设备标识 (UDI)

**端点**：`https://api.fda.gov/device/udi.json`

**用途**：访问包含设备标识信息的全球唯一设备标识数据库 (GUDID)。

**数据来源**：GUDID

**关键字段**：
- `identifiers.id` - 设备标识符 (DI)
- `identifiers.issuing_agency` - 签发机构（GS1、HIBCC、ICCBBA）
- `identifiers.type` - 主要或包 DI
- `brand_name` - 品牌名称
- `version_model_number` - 版本/型号
- `catalog_number` - 目录号
- `company_name` - 设备公司
- `device_count_in_base_package` - 基础包装中的数量
- `device_description` - 描述
- `is_rx` - 处方设备（真/假）
- `is_otc` - 非处方设备（真/假）
- `is_combination_product` - 组合产品（真/假）
- `is_kit` - 套件（真/假）
- `is_labeled_no_nrl` - 无乳胶标签
- `has_lot_or_batch_number` - 使用批号
- `has_serial_number` - 使用序列号
- `has_manufacturing_date` - 有生产日期
- `has_expiration_date` - 有过期日期
- `mri_safety` - MRI 安全状态
- `gmdn_terms` - 全球医疗器械命名术语
- `product_codes` - FDA 产品代码
- `storage` - 存储要求
- `customer_contacts` - 联系信息

**常见用例**：
- 设备识别和验证
- 供应链追踪
- 不良事件报告
- 库存管理
- 采购

**查询示例**：
```python
# Look up device by UDI
params = {
    "api_key": api_key,
    "search": "identifiers.id:00884838003019",
    "limit": 1
}

response = requests.get("https://api.fda.gov/device/udi.json", params=params)
```

```python
# Find prescription devices by brand name
params = {
    "api_key": api_key,
    "search": "brand_name:*insulin+pump*+AND+is_rx:true",
    "limit": 10
}
```

```python
# Search for MRI safe devices
params = {
    "api_key": api_key,
    "search": 'mri_safety:"MR Safe"',
    "limit": 50
}
```

### 9. COVID-19 血清学检测评估

**端点**：`https://api.fda.gov/device/covid19serology.json`

**目的**：获取 FDA 对 COVID-19 抗体测试的独立评估。

**数据来源**：FDA COVID-19 血清学测试性能

**关键字段**：
- `manufacturer` - 测试制造商
- `device` - 设备/测试名称
- `authorization_status` - EUA 状态
- `control_panel` - 用于评估的控制面板
- `sample_sensitivity_report_one` - 敏感度数据（第一份报告）
- `sample_specificity_report_one` - 特异性数据（第一份报告）
- `sample_sensitivity_report_two` - 敏感度数据（第二份报告）
- `sample_specificity_report_two` - 特异性数据（第二份报告）

**常见用例**：
- 测试性能比较
- 诊断准确性评估
- 采购决策支持
- 质量保证

**查询示例**：
```python
# Find tests by manufacturer
params = {
    "api_key": api_key,
    "search": "manufacturer:Abbott",
    "limit": 10
}

response = requests.get("https://api.fda.gov/device/covid19serology.json", params=params)
```

```python
# Get all tests with EUA
params = {
    "api_key": api_key,
    "search": "authorization_status:*EUA*",
    "limit": 100
}
```

## 集成技巧

### 全面的设备搜索

```python
def search_device_across_databases(device_name, api_key):
    """
    Search for a device across multiple FDA databases.

    Args:
        device_name: Name or partial name of device
        api_key: FDA API key

    Returns:
        Dictionary with results from each database
    """
    results = {}

    # Search adverse events
    events_url = "https://api.fda.gov/device/event.json"
    events_params = {
        "api_key": api_key,
        "search": f"device.brand_name:*{device_name}*",
        "limit": 10
    }
    results["adverse_events"] = requests.get(events_url, params=events_params).json()

    # Search 510(k) clearances
    fiveten_url = "https://api.fda.gov/device/510k.json"
    fiveten_params = {
        "api_key": api_key,
        "search": f"device_name:*{device_name}*",
        "limit": 10
    }
    results["510k_clearances"] = requests.get(fiveten_url, params=fiveten_params).json()

    # Search recalls
    recall_url = "https://api.fda.gov/device/enforcement.json"
    recall_params = {
        "api_key": api_key,
        "search": f"product_description:*{device_name}*",
        "limit": 10
    }
    results["recalls"] = requests.get(recall_url, params=recall_params).json()

    # Search UDI
    udi_url = "https://api.fda.gov/device/udi.json"
    udi_params = {
        "api_key": api_key,
        "search": f"brand_name:*{device_name}*",
        "limit": 10
    }
    results["udi"] = requests.get(udi_url, params=udi_params).json()

    return results
```

### 产品代码查找

```python
def get_device_classification(product_code, api_key):
    """
    Get detailed classification information for a device product code.

    Args:
        product_code: Three-letter FDA product code
        api_key: FDA API key

    Returns:
        Classification details dictionary
    """
    url = "https://api.fda.gov/device/classification.json"
    params = {
        "api_key": api_key,
        "search": f"product_code:{product_code}",
        "limit": 1
    }

    response = requests.get(url, params=params)
    data = response.json()

    if "results" in data and len(data["results"]) > 0:
        classification = data["results"][0]
        return {
            "product_code": classification.get("product_code"),
            "device_name": classification.get("device_name"),
            "device_class": classification.get("device_class"),
            "regulation_number": classification.get("regulation_number"),
            "medical_specialty": classification.get("medical_specialty_description"),
            "gmp_exempt": classification.get("gmp_exempt_flag") == "Y",
            "implant": classification.get("implant_flag") == "Y",
            "life_sustaining": classification.get("life_sustain_support_flag") == "Y"
        }
    return None
```

## 最佳实践

1. **使用产品代码** - 跨设备数据库搜索的最有效方法
2. **检查多个数据库** - 设备信息分布在多个端点
3. **处理大型结果集** - 设备数据库可能非常大；使用分页
4. **验证设备标识符** - 确保 UDI、510(k) 编号和 PMA 编号格式正确
5. **按设备类别过滤** - 在相关时按风险分类缩小搜索范围
6. **使用精确的品牌名称** - 通配符有效，但精确匹配更可靠
7. **考虑日期范围** - 设备数据积累数十年；适当时按日期过滤
8. **交叉参考数据** - 将不良事件与召回和注册联系起来以获取完整信息
9. **监控召回状态** - 召回状态从“正在进行”更改为“已完成”
10. **检查机构注册** - 机构必须每年注册一次；检查有效期

## 其他资源

- OpenFDA 设备 API 文档：https://open.fda.gov/apis/device/
- 设备分类数据库：https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfpcd/classification.cfm
- GUDID：https://accessgudid.nlm.nih.gov/
- API 基础知识：请参阅此参考目录中的 `api_basics.md`
- Python 示例：参见 `scripts/fda_device_query.py`