<!-- 此文件由机器翻译自 api_reference.md -->

# 开放目标平台 API 参考

## API 端点

```
https://api.platform.opentargets.org/api/v4/graphql
```

带有文档的交互式 GraphQL 游乐场：
<<<代码块_1>>>

## 访问方法

开放目标平台提供多种访问方式：

1. **GraphQL API** - 最适合单个实体查询和灵活的数据检索
2. **Web 界面** - 互动平台 https://platform.opentargets.org
3. **数据下载** - FTP 位于 https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/
4. **Google BigQuery** - 用于大规模系统查询

## 身份验证

GraphQL API 不需要身份验证。所有数据均可免费访问。

## 速率限制

对于涉及多个目标或疾病的系统查询，请使用数据集下载或 BigQuery，而不是重复的 API 调用。该 API 针对单实体和探索性查询进行了优化。

## GraphQL 查询结构

GraphQL 查询包括：
1.带有可选变量的查询操作
2.字段选择（只请求需要的字段）
3.嵌套实体遍历

### 基本 Python 示例

<<<代码块_2>>>

## 可用的查询端点

### /目标
检索基因注释、易处理性评估和疾病关联。

**常用字段：**
- `id` - Ensembl 基因 ID
- `approvedSymbol` - HGNC 基因符号
- `approvedName` - 完整基因名称
- `biotype` - 基因类型（蛋白质编码等）
- `tractability` - 成药性评估
- `safetyLiabilities` - 安全信息
- `expressions` - 基线表达数据
- `knownDrugs` - 批准/临床药物
- `associatedDiseases` - 疾病与证据的关联

### /疾病
检索疾病/表型数据、已知药物和临床信息。

**常用字段：**
- `id` - EFO 疾病标识符
- `name` - 疾病名称
- `description` - 疾病描述
- `therapeuticAreas` - 高级疾病类别
- `synonyms` - 备用名称
- `knownDrugs` - 适用于疾病的药物
- `associatedTargets` - 目标与证据的关联

### /药物
检索化合物详细信息、作用机制和药物警戒数据。

**常用字段：**
- `id` - ChEMBL 标识符
- `name` - 药物名称
- `drugType` - 小分子、抗体等。
- `maximumClinicalTrialPhase` - 开发阶段
- `indications` - 疾病适应症
- `mechanismsOfAction` - 目标机制
- `adverseEvents` - 药物警戒数据

### /搜索
搜索所有实体（目标、疾病、药物）。

**参数：**
- `queryString` - 搜索词
- `entityNames` - 按实体类型过滤
- `page` - 分页

### /关联疾病间接
检索目标疾病关联，包括来自本体论中疾病后代的间接证据。

**关键字段：**
- `rows` - 将记录与分数关联
- `aggregations` - 聚合统计信息

## 查询示例

### 查询1：获取与疾病关联的目标信息

<<<代码块_3>>>

### 查询 2：搜索疾病

<<<代码块_4>>>

### 查询 3：获取目标疾病对的证据

<<<代码块_5>>>

### 查询 4：获取治疗某种疾病的已知药物

<<<代码块_6>>>

## 错误处理

即使出现错误，GraphQL 也会返回状态代码 200。检查响应结构：

```python
if 'errors' in response_data:
    print(f"GraphQL errors: {response_data['errors']}")
else:
    print(f"Data: {response_data['data']}")
```

## 最佳实践

1. **仅请求所需字段** - 最大限度地减少数据传输并缩短响应时间
2. **使用变量** - 使查询可重用且更安全
3. **处理分页** - 大多数列表字段支持使用 `page: {size: N, index: M}` 分页
4. **探索模式** - 使用 GraphQL 浏览器发现可用字段
5. **批量相关查询** - 如果可能，将多个实体获取合并在一个查询中
6. **缓存结果** - 将频繁访问的数据存储在本地以减少 API 调用
7. **使用 BigQuery 进行批量** - 切换到 BigQuery/下载进行系统分析

## 数据许可

所有开放目标平台数据均可免费获取。在研究或商业产品中使用数据时，请引用最新出版物：

奥乔亚，D.等人。 (2025) 开放目标平台：促进药物发现中治疗假设的建立。核酸研究，53(D1)：D1467-D1477。