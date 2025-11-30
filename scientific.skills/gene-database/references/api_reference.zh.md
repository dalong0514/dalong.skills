<!-- 此文件由机器翻译自 api_reference.md -->

# NCBI 基因 API 参考

本文档提供了用于以编程方式访问 NCBI 基因数据库的详细 API 文档。

## 目录

1. [电子实用程序 API](#e-utilities-api)
2. [NCBI 数据集 API](#ncbi-datasets-api)
3. [身份验证和速率限制](#authentication-and-rate-limits)
4. [错误处理](#error-handling)

---

## 电子实用程序 API

E-utilities（Entrez 编程实用程序）为 NCBI 的 Entrez 数据库提供稳定的接口。

### 基本网址

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/
```

### 常用参数

- `db` - 数据库名称（使用 `gene` 作为基因数据库）
- `api_key` - 用于更高速率限制的 API 密钥
- `retmode` - 返回格式（json、xml、文本）
- `retmax` - 返回的最大记录数

### ESearch - 搜索数据库

搜索与文本查询匹配的基因。

**端点：** `esearch.fcgi`

**参数：**
- `db=gene`（必需）- 要搜索的数据库
- `term`（必需）- 搜索查询
- `retmax` - 最大结果数（默认值：20）
- `retmode` - json 或 xml（默认值：xml）
- `usehistory=y` - 将大型结果集的结果存储在历史服务器上

**查询语法：**
- 基因符号：`BRCA1[gene]` 或 `BRCA1[gene name]`
- 生物体：`human[organism]` 或 `9606[taxid]`
- 组合术语：`BRCA1[gene] AND human[organism]`
- 疾病：`muscular dystrophy[disease]`
- 染色体：`17q21[chromosome]`
- GO 术语：`GO:0006915[biological process]`

**请求示例：**

<<<代码块_1>>>

**响应格式（JSON）：**

<<<代码块_2>>>

### ESummary - 文档摘要

检索基因 ID 的文档摘要。

**端点：** `esummary.fcgi`

**参数：**
- `db=gene`（必需）- 数据库
- `id`（必需） - 以逗号分隔的基因 ID（最多 500 个）
- `retmode` - json 或 xml（默认值：xml）

**请求示例：**

<<<代码块_3>>>

**响应格式（JSON）：**

<<<代码块_4>>>

### EFetch - 全记录

获取各种格式的详细基因记录。

**端点：** `efetch.fcgi`

**参数：**
- `db=gene`（必需）- 数据库
- `id`（必需）- 逗号分隔的基因 ID
- `retmode` - xml、文本、asn.1（默认值：xml）
- `rettype` - 基因表、文档和

**请求示例：**

<<<代码块_5>>>

**XML 响应：** 包含详细的基因信息，包括：
- 基因命名法
- 序列位置
- 转录变体
- 蛋白质产品
- 基因本体注释
- 交叉引用
- 出版物

### ELink - 相关记录

在 Gene 或其他数据库中查找相关记录。

**端点：** `elink.fcgi`

**参数：**
- `dbfrom=gene`（必需）- 源数据库
- `db`（必需） - 目标数据库（基因、nuccore、蛋白质、pubmed 等）
- `id`（必需）- 基因 ID

**请求示例：**

<<<代码块_6>>>

### EInfo - 数据库信息

获取有关基因数据库的信息。

**端点：** `einfo.fcgi`

**参数：**
- `db=gene` - 要查询的数据库

**请求示例：**

```bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=gene&retmode=json"
```

---

## NCBI 数据集 API

数据集 API 通过元数据和序列提供对基因数据的简化访问。

### 基本网址

```
https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene
```

### 身份验证

在请求标头中包含 API 密钥：

```
api-key: YOUR_API_KEY
```

### 通过 ID 获取基因

通过基因 ID 检索基因数据。

**端点：** `GET /gene/id/{gene_id}`

**请求示例：**

```bash
curl "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/id/672"
```

**响应格式（JSON）：**

```json
{
  "genes": [
    {
      "gene": {
        "gene_id": "672",
        "symbol": "BRCA1",
        "description": "BRCA1 DNA repair associated",
        "tax_name": "Homo sapiens",
        "taxid": 9606,
        "chromosomes": ["17"],
        "type": "protein-coding",
        "synonyms": ["BRCC1", "FANCS", "PNCA4", "RNF53"],
        "nomenclature_authority": {
          "authority": "HGNC",
          "identifier": "HGNC:1100"
        },
        "genomic_ranges": [
          {
            "accession_version": "NC_000017.11",
            "range": [
              {
                "begin": 43044295,
                "end": 43170245,
                "orientation": "minus"
              }
            ]
          }
        ],
        "transcripts": [
          {
            "accession_version": "NM_007294.4",
            "length": 7207
          }
        ]
      }
    }
  ]
}
```

### 通过符号获取基因

按符号和生物体检索基因数据。

**端点：** `GET /gene/symbol/{symbol}/taxon/{taxon}`

**参数：**
- `{symbol}` - 基因符号（例如，BRCA1）
- `{taxon}` - 分类单元 ID（例如，人类为 9606）

**请求示例：**

```bash
curl "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/symbol/BRCA1/taxon/9606"
```

### 获取多个基因

检索多个基因的数据。

**端点：** `POST /gene/id`

**请求正文：**

```json
{
  "gene_ids": ["672", "7157", "5594"]
}
```

**请求示例：**

```bash
curl -X POST "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/gene/id" \
  -H "Content-Type: application/json" \
  -d '{"gene_ids": ["672", "7157", "5594"]}'
```

---

## 身份验证和速率限制

### 获取 API 密钥

1. 在 https://www.ncbi.nlm.nih.gov/account/ 创建 NCBI 帐户
2. 导航至设置 → API 密钥管理
3. 生成新的API密钥
4. 在请求中包含密钥

### 速率限制

**电子公用事业：**
- 没有 API 密钥：3 个请求/秒
- 使用 API 密钥：10 个请求/秒

**数据集API：**
- 没有 API 密钥：5 个请求/秒
- 使用 API 密钥：10 个请求/秒

### 使用指南

1. **在请求中包含电子邮件：** 将 `&email=your@email.com` 添加到电子公用事业请求中
2. **实施速率限制：** 使用请求之间的延迟
3. **使用 POST 进行大型查询：** 当使用许多 ID 时
4. **缓存结果：** 将经常访问的数据存储在本地
5. **优雅地处理错误：** 使用指数退避实现重试逻辑

---

## 错误处理

### HTTP 状态代码

- `200 OK` - 请求成功
- `400 Bad Request` - 参数无效或查询格式错误
- `404 Not Found` - 未找到基因 ID 或符号
- `429 Too Many Requests` - 超出速率限制
- `500 Internal Server Error` - 服务器错误（退避重试）

### 电子实用程序错误消息

电子实用程序在响应正文中返回错误：

**XML 格式：**
```xml
<ERROR>Empty id list - nothing to do</ERROR>
```

**JSON 格式：**
```json
{
  "error": "Invalid db name"
}
```

### 常见错误

1. **空结果集**
   - 原因：未找到基因符号或ID
   - 解决方案：验证拼写，检查生物过滤器

2. **超出速率限制**
   - 原因：请求太多
   - 解决方案：添加延迟，使用 API 密钥

3. **无效的查询语法**
   - 原因：搜索词格式错误
   - 解决方案：使用正确的字段标签（例如，`[gene]`、`[organism]`）

4. **超时**
   - 原因：结果集较大或连接速度较慢
   - 解决方案：使用历史记录服务器，减少结果大小

### 重试策略

对失败的请求实施指数退避：

```python
import time

def retry_request(func, max_attempts=3):
    for attempt in range(max_attempts):
        try:
            return func()
        except Exception as e:
            if attempt < max_attempts - 1:
                wait_time = 2 ** attempt  # 1s, 2s, 4s
                time.sleep(wait_time)
            else:
                raise
```

---

## 常见分类单元 ID

|有机体 |学名|分类单元 ID |
|----------|----------------|----------|
|人类 |智人 | 9606 | 9606
|鼠标|小家鼠 | 10090 | 10090
|老鼠 |褐家鼠 | 10116 |
|斑马鱼 |斑马鱼 | 7955 | 7955
|果蝇|果蝇 | 7227 | 7227
|线虫 |秀丽隐杆线虫 | 6239 | 6239
|酵母|酿酒酵母| 4932 |
|拟南芥|拟南芥| 3702 | 3702
|大肠杆菌 |大肠杆菌 | 562 | 562

---

## 其他资源

- **电子实用程序文档：** https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **数据集 API 文档：** https://www.ncbi.nlm.nih.gov/datasets/docs/v2/
- **基因数据库帮助：** https://www.ncbi.nlm.nih.gov/gene/
- **API 密钥注册：** https://www.ncbi.nlm.nih.gov/account/