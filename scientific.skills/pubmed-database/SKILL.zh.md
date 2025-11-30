<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：pubmed 数据库
描述：“直接 REST API 访问 PubMed。高级布尔/MeSH 查询、电子实用程序 API、批处理、引文管理。对于 Python 工作流程，更喜欢 biopython (Bio.Entrez)。使用它来直接 HTTP/REST 工作或自定义 API 实现。”
---

# PubMed 数据库

## 概述

PubMed 是美国国家医学图书馆的综合数据库，提供免费访问 MEDLINE 和生命科学文献。使用布尔运算符、MeSH 术语和字段标签构建高级查询，通过 E-utilities API 以编程方式访问数据，以进行系统评价和文献分析。

## 何时使用此技能

该技能应该在以下情况下使用：
- 搜索生物医学或生命科学研究文章
- 使用布尔运算符、字段标签或 MeSH 术语构建复杂的搜索查询
- 进行系统的文献综述或荟萃分析
- 通过电子实用程序 API 以编程方式访问 PubMed 数据
- 按特定标准查找文章（作者、期刊、出版日期、文章类型）
- 检索引文信息、摘要或全文文章
- 使用 PMID (PubMed ID) 或 DOI
- 创建用于文献监控或数据提取的自动化工作流程

## 核心能力

### 1. 高级搜索查询构建

使用布尔运算符、字段标签和专用语法构建复杂的 PubMed 查询。

**基本搜索策略**：
- 将概念与布尔运算符（AND、OR、NOT）结合起来
- 使用字段标签将搜索限制为特定记录部分
- 使用双引号进行短语搜索以实现精确匹配
- 对术语变体应用通配符
- 使用邻近搜索指定距离内的术语

**查询示例**：
```
# Recent systematic reviews on diabetes treatment
diabetes mellitus[mh] AND treatment[tiab] AND systematic review[pt] AND 2023:2024[dp]

# Clinical trials comparing two drugs
(metformin[nm] OR insulin[nm]) AND diabetes mellitus, type 2[mh] AND randomized controlled trial[pt]

# Author-specific research
smith ja[au] AND cancer[tiab] AND 2023[dp] AND english[la]
```

**何时查阅 search_syntax.md**：
- 需要可用字段标签的完整列表
- 需要搜索运算符的详细解释
- 构建复杂的邻近搜索
- 了解自动术语映射行为
- 需要日期范围、通配符或特殊字符的特定语法

字段标签的 Grep 模式：`\[au\]|\[ti\]|\[ab\]|\[mh\]|\[pt\]|\[dp\]`

### 2. MeSH 术语和受控词汇

使用医学主题词 (MeSH) 在生物医学文献中进行精确、一致的搜索。

**MeSH 搜索**：
- [mh] 标签搜索 MeSH 术语并自动包含较窄术语
- [majr] 标签限制为以主题为主要焦点的文章
- 将 MeSH 术语与副标题相结合以提高特异性（例如，糖尿病/治疗[mh]）

**常见 MeSH 副标题**：
- /diagnosis - 诊断方法
- /药物治疗 - 药物治疗
- /流行病学 - 疾病模式和患病率
- /病因学 - 疾病原因
- /预防与控制 - 预防措施
- /therapy - 治疗方法

**示例**：
<<<代码块_1>>>

### 3. 文章类型和出版物过滤

按出版物类型、日期、文本可用性和其他属性过滤结果。

**发布类型**（使用 [pt] 字段标签）：
- 临床试验
- 荟萃分析
- 随机对照试验
- 回顾
- 系统审查
- 病例报告
- 指南

**日期过滤**：
- 单年：`2024[dp]`
- 日期范围：`2020:2024[dp]`
- 具体日期：`2024/03/15[dp]`

**文本可用性**：
- 免费全文：添加`AND free full text[sb]`进行查询
- 有摘要：添加`AND hasabstract[text]`进行查询

**示例**：
<<<代码块_2>>>

### 4. 通过电子实用程序 API 进行编程访问

使用 NCBI 电子实用程序 REST API 以编程方式访问 PubMed 数据，以实现自动化和批量操作。

**核心 API 端点**：
1. **ESearch** - 搜索数据库并检索 PMID
2. **EFetch** - 下载各种格式的完整记录
3. **ESummary** - 获取文档摘要
4. **EPost** - 上传UID进行批量处理
5. **ELink** - 查找相关文章和链接数据

**基本工作流程**：
<<<代码块_3>>>

**速率限制**：
- 没有 API 密钥：3 个请求/秒
- 使用 API 密钥：10 个请求/秒
- 始终包含 User-Agent 标头

**最佳实践**：
- 对大型结果集使用历史服务器 (usehistory=y)
- 通过EPost实现多个UID的批量操作
- 在本地缓存结果以最大程度地减少冗余调用
- 遵守速率限制以避免服务中断

**何时查阅 api_reference.md**：
- 需要详细的端点文档
- 需要每个电子实用程序的参数规范
- 构建批量操作或历史服务器工作流程
- 了解响应格式（XML、JSON、文本）
- 解决 API 错误或速率限制问题

API 端点的 Grep 模式：`esearch|efetch|esummary|epost|elink|einfo`

### 5. 引文匹配和文章检索

使用部分引用信息或特定标识符查找文章。

**按标识符**：
<<<代码块_4>>>

**引文匹配**（通过 ECitMatch API）：
使用期刊名称、年份、卷数、页数和作者来查找 PMID：
<<<代码块_5>>>

**按作者和元数据**：
<<<代码块_6>>>

### 6.系统文献综述

进行全面的文献检索以进行系统评价和荟萃分析。

**PICO 框架**（人口、干预、比较、结果）：
系统地构建临床研究问题：
```
# Example: Diabetes treatment effectiveness
# P: diabetes mellitus, type 2[mh]
# I: metformin[nm]
# C: lifestyle modification[tiab]
# O: glycemic control[tiab]

diabetes mellitus, type 2[mh] AND
(metformin[nm] OR lifestyle modification[tiab]) AND
glycemic control[tiab] AND
randomized controlled trial[pt]
```

**综合搜索策略**：
```
# Include multiple synonyms and MeSH terms
(disease name[tiab] OR disease name[mh] OR synonym[tiab]) AND
(treatment[tiab] OR therapy[tiab] OR intervention[tiab]) AND
(systematic review[pt] OR meta-analysis[pt] OR randomized controlled trial[pt]) AND
2020:2024[dp] AND
english[la]
```

**搜索细化**：
1. 广泛着手，审查结果
2. 使用字段标签添加特异性
3.应用日期和出版物类型过滤器
4. 使用高级搜索查看查询翻译
5. 合并复杂查询的搜索历史

**何时查阅 common_queries.md**：
- 需要特定疾病类型或研究领域的示例查询
- 需要不同研究设计的模板
- 寻找特定人群的查询模式（儿科、老年人等）
- 构建特定于方法的搜索
- 需要质量过滤器或最佳实践模式

查询示例的 Grep 模式：`diabetes|cancer|cardiovascular|clinical trial|systematic review`

### 7. 搜索历史和保存的搜索

使用 PubMed 的检索历史和 My NCBI 功能实现高效的研究工作流程。

**搜索历史**（通过高级搜索）：
- 维护最多 100 条搜索
- 8 小时不活动后过期
- 使用 # 个参考文献合并之前的搜索
- 执行前预览结果计数

**示例**：
```
#1: diabetes mellitus[mh]
#2: cardiovascular diseases[mh]
#3: #1 AND #2 AND risk factors[tiab]
```

**我的 NCBI 功能**：
- 无限期保存搜索
- 为新的匹配文章设置电子邮件提醒
- 创建已保存文章的集合
- 按项目或主题组织研究

**RSS 提要**：
为任何搜索创建 RSS 源，以监控您感兴趣领域的新出版物。

### 8. 相关文章和引文发现

查找相关研究并探索引用网络。

**类似文章功能**：
每篇 PubMed 文章都包含基于以下内容预先计算的相关文章：
- 标题和摘要的相似性
- MeSH 术语重叠
- 加权算法匹配

**相关数据ELink**：
```
# Find related articles programmatically
elink.fcgi?dbfrom=pubmed&db=pubmed&id=PMID&cmd=neighbor
```

**引文链接**：
- 链接到出版商的全文
- PubMed Central 免费文章的链接
- 与相关 NCBI 数据库（GenBank、ClinicalTrials.gov 等）的连接

### 9. 导出和引文管理

以各种格式导出搜索结果以进行引文管理和进一步分析。

**导出格式**：
- 用于参考管理器的 .nbib 文件（Zotero、Mendeley、EndNote）
- AMA、MLA、APA、NLM 引文样式
- 用于数据分析的CSV
- 用于编程处理的 XML

**剪贴板和收藏夹**：
- 剪贴板：临时存储最多 500 个项目（8 小时有效期）
- 馆藏：通过我的 NCBI 帐户永久存储

**通过API批量导出**：
```python
# Export citations in MEDLINE format
efetch.fcgi?db=pubmed&id=PMID1,PMID2&rettype=medline&retmode=text
```

## 使用参考文件

该技能包括`references/`目录中的三个综合参考文件：

### 参考文献/api_reference.md
完整的电子公用事业 API 文档，包括所有九个端点、参数、响应格式和最佳实践。咨询时间：
- 实施程序化的 PubMed 访问
- 构建API请求
- 了解速率限制和身份验证
- 通过历史服务器处理大型数据集
- 排查 API 错误

### 参考文献/search_syntax.md
PubMed 搜索语法的详细指南，包括字段标签、布尔运算符、通配符和特殊字符。咨询时间：
- 构建复杂的搜索查询
- 了解自动术语映射
- 使用高级搜索功能（邻近度、通配符）
- 应用过滤器和限制
- 排除意外搜索结果的故障

### 参考文献/common_queries.md
广泛收集各种研究场景、疾病类型和方法的示例查询。咨询时间：
- 开始新的文献检索
- 需要特定研究领域的模板
- 寻找最佳实践查询模式
- 进行系统审查
- 搜索特定的研究设计或人群

**参考加载策略**：
根据特定任务的需要，将参考文件加载到上下文中。对于简短查询或基本搜索，此 SKILL.md 中的信息可能就足够了。对于复杂的操作，请查阅相应的参考文件。

## 常见工作流程

### 工作流程 1：基础文献检索

1. 识别关键概念和同义词
2. 使用布尔运算符和字段标签构建查询
3. 检查初始结果并优化查询
4.应用过滤器（日期、文章类型、语言）
5. 导出结果进行分析

### 工作流程 2：系统评论搜索

1.使用PICO框架定义研究问题
2. 识别所有相关的 MeSH 术语和同义词
3. 构建综合搜索策略
4. 检索多个数据库（包括PubMed）
5. 文献检索策略及日期
6. 导出结果进行筛选和审核

### 工作流程 3：程序化数据提取

1.设计搜索查询并在Web界面中进行测试
2.使用ESearch API实现搜索
3. 使用历史服务器处理大型结果集
4. 使用 EFetch 检索详细记录
5. 解析 XML/JSON 响应
6. 使用缓存将数据存储在本地
7. 实施速率限制和错误处理

### 工作流程 4：引文发现

1. 从已知的相关文章开始
2. 使用类似文章查找相关工作
3. 检查引用的文章（如果有）
4. 从相关文章中探索 MeSH 术语
5. 根据发现构建新的搜索
6.使用ELink查找相关数据库条目

### 工作流程 5：持续文献监测

1. 构建综合搜索查询
2. 测试和细化查询的精度
3. 将搜索保存到我的 NCBI 帐户
4. 为新比赛设置电子邮件提醒
5. 创建 RSS feed 用于 feed 阅读器监控
6.定期回顾新文章

## 提示和最佳实践

### 搜索策略
- 从广泛开始，然后使用字段标签和过滤器缩小范围
- 包括同义词和 MeSH 术语以实现全面覆盖
- 对确切的短语使用引号
- 检查高级搜索中的搜索详细信息以验证查询翻译
- 使用搜索历史合并多个搜索

### API 使用
- 获取 API 密钥以获得更高的速率限制（10 请求/秒 vs 3 请求/秒）
- 使用历史服务器获取超过 500 篇文章的结果集
- 实施指数退避以进行速率限制处理
- 在本地缓存结果以最大程度地减少冗余请求
- 始终包含描述性的 User-Agent 标头

### 质量过滤
- 更喜欢对综合证据进行系统评价和荟萃分析
- 使用出版物类型过滤器来查找特定的研究设计
- 按日期过滤最新研究
- 根据需要应用语言过滤器
- 使用免费的全文过滤器立即访问

### 引文管理
- 尽早并经常导出以避免丢失搜索结果
- 使用 .nbib 格式与大多数参考管理器兼容
- 为永久馆藏创建我的 NCBI 帐户
- 可重复性的文档搜索策略
- 使用集合按项目组织研究

## 限制和注意事项

### 数据库覆盖范围
- 主要是生物医学和生命科学文献
- 1975 年之前的文章通常缺乏摘要
- 2002 年以后的完整作者姓名
- 提供非英文摘要，但可能默认显示英文

### 搜索限制
- 最多显示 10,000 个结果
- 搜索历史记录在 8 小时不活动后过期
- 剪贴板最多可容纳 500 个项目，有效期为 8 小时
- 自动术语映射可能会产生意想不到的结果

### API注意事项
- 适用速率限制（3-10 个请求/秒）
- 大型查询可能会超时（使用历史服务器）
- 详细数据提取需要 XML 解析
- 建议用于生产用途的 API 密钥

### 访问限制
- PubMed 提供引文和摘要（并不总是全文）
- 全文访问取决于出版商、机构访问或开放获取状态
- LinkOut 的可用性因期刊和机构而异
- 部分内容需要订阅或付费

## 支持资源

- **PubMed 帮助**：https://pubmed.ncbi.nlm.nih.gov/help/
- **电子实用程序文档**：https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **NLM 服务台**：1-888-FIND-NLM (1-888-346-3656)
- **技术支持**：vog.hin.mln.ibcn@seitilitue
- **邮件列表**：utilities-announce@ncbi.nlm.nih.gov