<!-- 此文件由机器翻译自 other.md -->

# FDA 其他数据库 - 物质和 NSDE

该参考文献涵盖了 FDA 物质相关的端点和可通过 openFDA 访问的其他专用 API 端点。

## 概述

FDA 维护着精确到分子水平的物质级信息的附加数据库。这些数据库支持药物、生物制剂、设备、食品和化妆品的监管活动。

## 可用端点

### 1. 物质数据

**端点**：`https://api.fda.gov/other/substance.json`

**目的**：获取精确到分子水平的物质信息，供内部和外部使用。这包括有关 FDA 监管产品中使用的活性药物成分、赋形剂和其他物质的信息。

**数据来源**：FDA 全球物质注册系统 (GSRS)

**关键字段**：
- `uuid` - 唯一物质标识符 (UUID)
- `approvalID` - FDA 唯一成分标识符 (UNII)
- `approved` - 批准日期
- `substanceClass` - 物质类型（化学品、蛋白质、核酸、聚合物等）
- `names` - 物质名称数组
- `names.name` - 名称文本
- `names.type` - 名称类型（系统、品牌、通用等）
- `names.preferred` - 是否首选名称
- `codes` - 物质代码数组
- `codes.code` - 代码值
- `codes.codeSystem` - 代码系统（CAS、ECHA、EINECS 等）
- `codes.type` - 代码类型
- `relationships` - 物质关系数组
- `relationships.type` - 关系类型（活性部分、代谢物、杂质等）
- `relationships.relatedSubstance` - 相关物质参考
- `moieties` - 分子部分
- `properties` - 物理化学特性数组
- `properties.name` - 属性名称
- `properties.value` - 属性值
- `properties.propertyType` - 属性类型
- `structure` - 化学结构信息
- `structure.smiles` - SMILES 表示法
- `structure.inchi` - InChI 字符串
- `structure.inchiKey` - InChI 键
- `structure.formula` - 分子式
- `structure.molecularWeight` - 分子量
- `modifications` - 结构修改（对于蛋白质等）
- `protein` - 蛋白质特定信息
- `protein.subunits` - 蛋白质亚基
- `protein.sequenceType` - 序列类型
- `nucleicAcid` - 核酸信息
- `nucleicAcid.subunits` - 序列子单元
- `polymer` - 聚合物信息
- `mixture` - 混合组件
- `mixture.components` - 成分物质
- `tags` - 物质标签
- `references` - 文献参考

**物质类别**：
- **化学** - 具有明确化学结构的小分子
- **蛋白质** - 蛋白质和肽
- **核酸** - DNA、RNA、寡核苷酸
- **聚合物** - 聚合物质
- **结构多样** - 复杂混合物、植物药
- **混合物** - 定义的混合物
- **概念** - 抽象概念（例如，组）

**常见用例**：
- 活性成分鉴定
- 分子结构查找
- UNII代码解析
- 化学标识符映射（CAS 到 UNII 等）
- 物质关系分析
- 辅料鉴定
- 植物物质信息
- 蛋白质和生物学表征

**查询示例**：
```python
import requests

api_key = "YOUR_API_KEY"
url = "https://api.fda.gov/other/substance.json"

# Look up substance by UNII code
params = {
    "api_key": api_key,
    "search": "approvalID:R16CO5Y76E",  # Aspirin UNII
    "limit": 1
}

response = requests.get(url, params=params)
data = response.json()
```

<<<代码块_1>>>

<<<代码块_2>>>

<<<代码块_3>>>

<<<代码块_4>>>

<<<代码块_5>>>

### 2. NSDE（国家物质数据库条目）

**端点**：`https://api.fda.gov/other/nsde.json`

**目的**：从旧版国家药品法规 (NDC) 目录条目中访问历史物质数据。该端点提供历史药品列表中出现的物质信息。

**注**：该数据库主要用于历史参考。有关当前物质信息，请使用物质数据端点。

**关键字段**：
- `proprietary_name` - 产品专有名称
- `nonproprietary_name` - 非专有名称
- `dosage_form` - 剂型
- `route` - 给药途径
- `company_name` - 公司名称
- `substance_name` - 物质名称
- `active_numerator_strength` - 活性成分强度（分子）
- `active_ingred_unit` - 活性成分单位
- `pharm_classes` - 药理学课程
- `dea_schedule` - DEA 管制物质表

**常见用例**：
- 历史药物配方研究
- 遗留系统集成
- 历史物质名称映射
- 制药史研究
**查询示例**：
<<<代码块_6>>>

```python
# Find controlled substances by DEA schedule
params = {
    "api_key": api_key,
    "search": "dea_schedule:CII",
    "limit": 50
}
```

## 集成技巧

### UNII 到 CAS 映射

```python
def get_substance_identifiers(unii, api_key):
    """
    Get all identifiers for a substance given its UNII code.

    Args:
        unii: FDA Unique Ingredient Identifier
        api_key: FDA API key

    Returns:
        Dictionary with substance identifiers
    """
    import requests

    url = "https://api.fda.gov/other/substance.json"
    params = {
        "api_key": api_key,
        "search": f"approvalID:{unii}",
        "limit": 1
    }

    response = requests.get(url, params=params)
    data = response.json()

    if "results" not in data or len(data["results"]) == 0:
        return None

    substance = data["results"][0]

    identifiers = {
        "unii": substance.get("approvalID"),
        "uuid": substance.get("uuid"),
        "preferred_name": None,
        "cas_numbers": [],
        "other_codes": {}
    }

    # Extract names
    if "names" in substance:
        for name in substance["names"]:
            if name.get("preferred"):
                identifiers["preferred_name"] = name.get("name")
                break
        if not identifiers["preferred_name"] and len(substance["names"]) > 0:
            identifiers["preferred_name"] = substance["names"][0].get("name")

    # Extract codes
    if "codes" in substance:
        for code in substance["codes"]:
            code_system = code.get("codeSystem", "").upper()
            code_value = code.get("code")

            if "CAS" in code_system:
                identifiers["cas_numbers"].append(code_value)
            else:
                if code_system not in identifiers["other_codes"]:
                    identifiers["other_codes"][code_system] = []
                identifiers["other_codes"][code_system].append(code_value)

    return identifiers
```

### 化学结构查找

```python
def get_chemical_structure(substance_name, api_key):
    """
    Get chemical structure information for a substance.

    Args:
        substance_name: Name of the substance
        api_key: FDA API key

    Returns:
        Dictionary with structure information
    """
    import requests

    url = "https://api.fda.gov/other/substance.json"
    params = {
        "api_key": api_key,
        "search": f"names.name:{substance_name}",
        "limit": 1
    }

    response = requests.get(url, params=params)
    data = response.json()

    if "results" not in data or len(data["results"]) == 0:
        return None

    substance = data["results"][0]

    if "structure" not in substance:
        return None

    structure = substance["structure"]

    return {
        "smiles": structure.get("smiles"),
        "inchi": structure.get("inchi"),
        "inchi_key": structure.get("inchiKey"),
        "formula": structure.get("formula"),
        "molecular_weight": structure.get("molecularWeight"),
        "substance_class": substance.get("substanceClass")
    }
```

### 物质关系图

```python
def get_substance_relationships(unii, api_key):
    """
    Get all related substances (metabolites, active moieties, etc.).

    Args:
        unii: FDA Unique Ingredient Identifier
        api_key: FDA API key

    Returns:
        Dictionary organizing relationships by type
    """
    import requests

    url = "https://api.fda.gov/other/substance.json"
    params = {
        "api_key": api_key,
        "search": f"approvalID:{unii}",
        "limit": 1
    }

    response = requests.get(url, params=params)
    data = response.json()

    if "results" not in data or len(data["results"]) == 0:
        return None

    substance = data["results"][0]

    relationships = {}

    if "relationships" in substance:
        for rel in substance["relationships"]:
            rel_type = rel.get("type")
            if rel_type not in relationships:
                relationships[rel_type] = []

            related = {
                "uuid": rel.get("relatedSubstance", {}).get("uuid"),
                "unii": rel.get("relatedSubstance", {}).get("approvalID"),
                "name": rel.get("relatedSubstance", {}).get("refPname")
            }
            relationships[rel_type].append(related)

    return relationships
```

### 活性成分提取

```python
def find_active_ingredients_by_product(product_name, api_key):
    """
    Find active ingredients in a drug product.

    Args:
        product_name: Drug product name
        api_key: FDA API key

    Returns:
        List of active ingredient UNIIs and names
    """
    import requests

    # First search drug label database
    label_url = "https://api.fda.gov/drug/label.json"
    label_params = {
        "api_key": api_key,
        "search": f"openfda.brand_name:{product_name}",
        "limit": 1
    }

    response = requests.get(label_url, params=label_params)
    data = response.json()

    if "results" not in data or len(data["results"]) == 0:
        return None

    label = data["results"][0]

    # Extract UNIIs from openfda section
    active_ingredients = []

    if "openfda" in label:
        openfda = label["openfda"]

        # Get UNIIs
        unii_list = openfda.get("unii", [])
        generic_names = openfda.get("generic_name", [])

        for i, unii in enumerate(unii_list):
            ingredient = {"unii": unii}
            if i < len(generic_names):
                ingredient["name"] = generic_names[i]

            # Get additional substance info
            substance_info = get_substance_identifiers(unii, api_key)
            if substance_info:
                ingredient.update(substance_info)

            active_ingredients.append(ingredient)

    return active_ingredients
```

## 最佳实践

1. **使用 UNII 作为主要标识符** - FDA 数据库中最一致
2. **标识符系统之间的映射** - CAS、UNII、InChI 用于交叉引用的密钥
3. **处理物质变化** - 不同的盐形式、水合物具有不同的 UNII
4. **检查物质类别** - 不同类别有不同的数据结构
5. **验证化学结构** - 应验证 SMILES 和 InChI
6. **考虑物质关系** - 活性部分与盐形式很重要
7. **使用首选名称** - 比商品名称更一致
8. **缓存物质数据** - 物质信息不经常更改
9. **与其他终点的交叉引用** - 将物质与药物/产品联系起来
10. **处理混合物成分** - 复杂产品具有多种成分

## UNII系统

FDA 唯一成分标识符 (UNII) 系统提供：
- **唯一标识符** - 每种物质都有一个 UNII
- **物质特异性** - 不同形式（盐、水合物）获得不同的 UNII
- **全球认可** - 在国际上使用
- **稳定性** - UNII 一旦分配就不会改变
- **免费访问** - 无需许可

**UNII 格式**：10 个字符的字母数字代码（例如，`R16CO5Y76E`）

## 物质类别解释

### 化学
- 传统小分子药物
- 已定义分子结构
- 包括有机和无机化合物
- SMILES、InChI、分子式可用

### 蛋白质
- 多肽和蛋白质
- 可用的序列信息
- 可能有翻译后修饰
- 包括抗体、酶、激素

### 核酸
- DNA和RNA序列
- 寡核苷酸
- 反义、siRNA、mRNA
- 可用的序列数据

### 聚合物
- 合成和天然聚合物
- 结构重复单元
- 分子量分布
- 用作赋形剂和活性成分

### 结构多样化
- 复杂的天然产物
- 植物提取物
- 无单分子结构的材料
- 以来源和成分为特征

### 混合物
- 定义的物质组合
- 固定或可变的成分
- 每个组件均可追踪

## 其他资源

- FDA 物质注册系统：https://fdasis.nlm.nih.gov/srs/
- UNII 搜索：https://precision.fda.gov/uniisearch
- OpenFDA 其他 API：https://open.fda.gov/apis/other/
- API 基础知识：请参阅此参考目录中的 `api_basics.md`
- Python 示例：参见 `scripts/fda_substance_query.py`