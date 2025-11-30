<!-- 此文件由机器翻译自 api_reference.md -->

# ChEMBL Web 服务 API 参考

## 概述

ChEMBL 是由欧洲生物信息学研究所 (EBI) 维护的手动管理的具有药物样特性的生物活性分子数据库。它包含有关化合物、靶标、分析、生物活性数据和批准药物的信息。

ChEMBL 数据库包含：
- 超过 200 万条化合物记录
- 超过 140 万条化验记录
- 超过 1900 万个活动值
- 超过 13,000 个药物靶点的信息
- 超过 16,000 种已批准药物和临床候选药物的数据

## Python 客户端安装

```bash
pip install chembl_webresource_client
```

## 关键资源和端点

ChEMBL 提供对 30 多个专用端点的访问：

### 核心数据类型

- **分子** - 化合物结构、性质和同义词
- **目标** - 蛋白质和非蛋白质生物目标
- **活性** - 生物测定测量结果
- **测定** - 实验测定详细信息
- **药物** - 批准的药物信息
- **机制** - 药物作用机制数据
- **文档** - 文献来源和参考文献
- **cell_line** - 细胞系信息
- **组织** - 组织类型
- **蛋白质_类别** - 蛋白质分类
- **target_component** - 目标组件详细信息
- **compound_structural_alert** - 毒性结构警报

## 查询模式和过滤器

### 过滤运算符

该 API 支持 Django 风格的过滤器运算符：

- `__exact` - 完全匹配
- `__iexact` - 不区分大小写的精确匹配
- `__contains` - 包含子字符串
- `__icontains` - 不区分大小写包含
- `__startswith` - 以前缀开头
- `__endswith` - 以后缀结尾
- `__gt` - 大于
- `__gte` - 大于或等于
- `__lt` - 小于
- `__lte` - 小于或等于
- `__range` - 范围内的值
- `__in` - 列表中的值
- `__isnull` - 为空/不为空
- `__regex` - 正则表达式匹配
- `__search` - 全文搜索

### 过滤器查询示例

**分子量过滤：**
<<<代码块_1>>>

**名称模式匹配：**
<<<代码块_2>>>

**多个条件：**
<<<代码块_3>>>

## 化学结构搜索

### 子结构搜索
使用 SMILES 搜索包含特定子结构的化合物：

<<<代码块_4>>>

### 相似性搜索
查找与查询结构相似的化合物：

<<<代码块_5>>>

## 常见数据检索模式

### 通过 ChEMBL ID 获取分子
<<<代码块_6>>>

### 获取目标信息
```python
target = new_client.target.get('CHEMBL240')
```

### 获取活动数据
```python
activities = new_client.activity.filter(
    target_chembl_id='CHEMBL240',
    standard_type='IC50',
    standard_value__lte=100
)
```

### 获取药物信息
```python
drug = new_client.drug.get('CHEMBL1234')
```

## 响应格式

API支持多种响应格式：
- JSON（默认）
- XML
-YAML

## 缓存和性能

Python客户端自动在本地缓存结果：
- **默认缓存持续时间**：24 小时
- **缓存位置**：本地文件系统
- **惰性求值**：仅在访问数据时才执行查询

### 配置设置

```python
from chembl_webresource_client.settings import Settings

# Disable caching
Settings.Instance().CACHING = False

# Adjust cache expiration (in seconds)
Settings.Instance().CACHE_EXPIRE = 86400  # 24 hours

# Set timeout
Settings.Instance().TIMEOUT = 30

# Set retries
Settings.Instance().TOTAL_RETRIES = 3
```

## 分子特性

可用的常见分子特性：

- `mw_freebase` - 分子量
- `alogp` - 计算的 LogP
- `hba` - 氢键受体
- `hbd` - 氢键供体
- `psa` - 极性表面积
- `rtb` - 可旋转键
- `ro3_pass` - 符合 3 规则
- `num_ro5_violations` - Lipinski 5 次违规规则
- `cx_most_apka` - 酸性最强的 pKa
- `cx_most_bpka` - 最基本的 pKa
- `molecular_species` - 分子种类
- `full_mwt` - 全分子量

## 生物活性数据字段

主要生物活性领域：

- `standard_type` - 活动类型（IC50、Ki、Kd、EC50 等）
- `standard_value` - 数字活动值
- `standard_units` - 单位（nM、uM 等）
- `pchembl_value` - 标准化活动值（-对数刻度）
- `activity_comment` - 活动注释
- `data_validity_comment` - 数据有效性标志
- `potential_duplicate` - 重复标志

## 目标信息字段

目标数据包括：

- `target_chembl_id` - ChEMBL 目标标识符
- `pref_name` - 首选目标名称
- `target_type` - 类型（蛋白质、生物体等）
- `organism` - 目标生物体
- `tax_id` - NCBI 分类 ID
- `target_components` - 组件详细信息

## 高级查询示例

### 寻找激酶抑制剂
```python
# Get kinase targets
targets = new_client.target.filter(
    target_type='SINGLE PROTEIN',
    pref_name__icontains='kinase'
)

# Get activities for these targets
activities = new_client.activity.filter(
    target_chembl_id__in=[t['target_chembl_id'] for t in targets],
    standard_type='IC50',
    standard_value__lte=100
)
```

### 检索药物机制
```python
mechanisms = new_client.mechanism.filter(
    molecule_chembl_id='CHEMBL25'
)
```
### 获取复合生物活性
```python
activities = new_client.activity.filter(
    molecule_chembl_id='CHEMBL25',
    pchembl_value__isnull=False
)
```

## 图像生成

ChEMBL 可以生成分子结构的 SVG 图像：

```python
from chembl_webresource_client.new_client import new_client
image = new_client.image
svg = image.get('CHEMBL25')
```

## 分页

结果自动分页。迭代所有结果：

```python
activities = new_client.activity.filter(target_chembl_id='CHEMBL240')
for activity in activities:
    print(activity)
```

## 错误处理

常见错误：
- **404**：找不到资源
- **503**：服务暂时不可用
- **超时**：请求花费的时间太长

客户端根据 `TOTAL_RETRIES` 设置自动重试失败的请求。

## 速率限制

ChEMBL 具有公平使用政策：
- 尊重查询频率
- 使用缓存来最大程度地减少重复请求
- 考虑批量下载大型数据集

## 其他资源

- 官方API文档：https://www.ebi.ac.uk/chembl/api/data/docs
- Python 客户端 GitHub：https://github.com/chembl/chembl_webresource_client
- ChEMBL 接口文档：https://chembl.gitbook.io/chembl-interface-documentation/
- 示例笔记本：https://github.com/chembl/notebooks