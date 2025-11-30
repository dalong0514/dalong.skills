<!-- 此文件由机器翻译自 additional_apis.md -->

# 其他 USPTO API 参考

## 概述

除了专利检索、PEDS 和商标之外，USPTO 还提供用于引文、审查意见通知书、转让、诉讼和其他专利数据的专用 API。

## 1. 丰富的引文 API

### 概述

深入了解专利评估流程并引用 IP5（美国专利商标局、欧洲专利局、日本特许厅、韩国特许厅、中国国家知识产权局）和公众使用的参考文献。

**版本：** v3、v2、v1

**基本 URL：** 通过 USPTO 开放数据门户访问

### 目的

分析审查员在专利审查期间引用了哪些参考文献以及专利如何引用现有技术。

### 主要特点

- **前向引用** - 引用给定专利的专利
- **向后引用** - 专利引用的参考文献
- **审查员引用** - 审查员与申请人引用的参考文献
- **引用上下文** - 引用参考文献的方式和原因

### 用例

- 现有技术分析
- 专利格局分析
- 识别相关技术
- 根据引用评估专利强度

## 2. Office Action API

### 2.1 Office Action 文本检索 API

**版本：** v1

### 目的

检索专利申请的完整全文审查意见函件文件。

### 特点

- 审查意见全文
- 限制、拒绝、反对
- 审查员修正
- 搜索信息

### 使用示例

```python
# Retrieve office action text by application number
def get_office_action_text(app_number, api_key):
    """
    Fetch full text of office actions for an application.
    Note: Integrate with PEDS to identify which office actions exist.
    """
    # API implementation
    pass
```

### 2.2 审查意见通知书 API

**版本：** v2、测试版 v1

### 目的

提供从审查意见通知书中提取的专利引文数据，显示审查员在审查过程中使用的参考文献。

### 关键数据

- 专利和非专利文献引用
- 引文上下文（拒绝、信息等）
- 审查员搜索策略
- 检方研究数据集

### 2.3 审查意见驳回 API

**版本：** v2、测试版 v1

### 目的

详细说明拒绝原因和审查结果以及截至 2025 年 3 月的批量拒绝数据。

### 拒绝类型

- **35 U.S.C. § 102** - 预期（缺乏新颖性）
- **35 U.S.C. § 103** - 显而易见性
- **35 U.S.C. § 112** - 启用、书面描述、不确定性
- **35 U.S.C. § 101** - 主题资格

### 用例

- 分析常见的拒绝原因
- 识别有问题的索赔语言
- 根据历史数据准备回复
- 拒绝模式的投资组合分析

### 2.4 Office Action 每周 Zips API

**版本：** v1

### 目的

提供按每周发布时间表组织的办公行动文档全文的批量下载。

### 特点

- 每周档案下载
- 完整的办公行动文本
- 批量访问以进行大规模分析

## 3.专利转让检索API

### 概述

**版本：** v1.4

访问 USPTO 专利转让数据库以获取所有权记录和转让。

**基本网址：** `https://assignment-api.uspto.gov/patent/`

### 目的

跟踪专利所有权、转让、担保权益和公司交易。

### 搜索方法

#### 按专利号

<<<代码块_1>>>

#### 按申请号

<<<代码块_2>>>

#### 按受让人姓名

<<<代码块_3>>>

### 响应格式

返回带有类似于商标转让的转让记录的 XML：

- 卷轴/框架编号
- 输送类型
- 日期（执行和记录）
- 转让人和受让人
- 受影响的专利/申请

### 常见用途

<<<代码块_4>>>

## 4. PTAB API（专利审判和上诉委员会）

### 概述

**版本：** v2

获取专利审判和上诉委员会诉讼数据。

### 目的

检索有关以下内容的信息：
- 多方复审（IPR）
- 授权后审查（PGR）
- 涵盖业务方法（CBM）审查
- 单方面上诉

### 可用数据

- 请愿信息
- 审判决定
- 最终书面决定
- 申请人和专利所有人信息
- 索赔受到质疑
- 试验结果

### 注意

目前正在迁移到新的开放数据门户。检查当前文档以获取访问详细信息。

## 5.专利诉讼案件API

### 概述

**版本：** v1

包含超过 74,623 条地区法院诉讼记录，涵盖专利诉讼数据。

### 目的

查看联邦地区法院的专利侵权案件。

### 关键数据

- 案件编号和提交日期
- 已主张的专利
- 当事人（原告和被告）
- 场地
- 案件结果

### 用例

- 诉讼风险分析
- 识别经常提起诉讼的专利
- 跟踪诉讼趋势
- 分析场地偏好
- 评估专利执法模式

## 6. Cancer Moonshot 专利数据集 API

### 概述

**版本：** v1.0.1

癌症相关专利发现的专业数据集。
### 目的

搜索和下载与癌症研究、治疗和诊断相关的专利。

### 特点

- 策划癌症相关专利
- 批量数据下载
- 按癌症类型分类
- 治疗方式分类

### 用例

- 癌症研究现有技术
- 技术格局分析
- 确定研究趋势
- 许可机会

## 7. OCE 专利审查状态/事件代码 API

### 概述

**版本：** v1

提供专利审查中使用的 USPTO 状态和事件代码的官方描述。

### 目的

解码 PEDS 和其他检查数据中的事务代码和状态代码。

### 提供的数据

- **状态代码** - 应用程序状态描述
- **事件代码** - 交易/事件描述
- **代码定义** - 官方含义

### 整合

与 PEDS 数据一起使用来解释交易代码：

<<<代码块_5>>>

## API 集成模式

### 组合工作流程示例

<<<代码块_6>>>

### 投资组合分析示例

```python
def analyze_company_portfolio(company_name, api_key):
    """
    Analyze a company's patent portfolio using multiple APIs.
    """
    # 1. Find all assigned patents
    assignments = find_company_patents(company_name, api_key)
    patent_numbers = extract_patent_numbers(assignments)

    # 2. Get details for each patent
    portfolio = []
    for patent_num in patent_numbers:
        patent_data = {
            'number': patent_num,
            'details': search_patent(patent_num, api_key),
            'citations': get_citations(patent_num, api_key),
            'litigation': get_litigation(patent_num, api_key)
        }
        portfolio.append(patent_data)

    # 3. Aggregate statistics
    stats = {
        'total_patents': len(portfolio),
        'cited_by_count': sum(len(p['citations']) for p in portfolio),
        'litigated_count': sum(1 for p in portfolio if p['litigation']),
        'technology_areas': aggregate_tech_areas(portfolio)
    }

    return {'portfolio': portfolio, 'statistics': stats}
```

## 最佳实践

1. **API密钥管理** - 使用环境变量，从不硬编码
2. **速率限制** - 对所有 API 实施指数退避
3. **缓存** - 缓存 API 响应以最大程度地减少冗余调用
4. **错误处理** - 优雅地处理 API 错误和丢失数据
5. **数据验证** - 在 API 调用之前验证输入格式
6. **组合 API** - 结合使用适当的 API 进行综合分析
7. **文档** - 跟踪 API 版本和更改

## API 密钥注册

所有 API 都需要在以下位置注册：
**https://account.uspto.gov/api-manager/**

单个 API 密钥适用于大多数 USPTO API。

## 资源

- **开发者门户**：https://developer.uspto.gov/
- **打开数据门户**：https://data.uspto.gov/
- **API目录**：https://developer.uspto.gov/api-catalog
- **Swagger 文档**：可用于单个 API