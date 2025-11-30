<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：基因数据库
描述：“通过电子实用程序/数据集 API 查询 NCBI 基因。按符号/ID 搜索，检索基因信息（RefSeqs、GO、位置、表型）、批量查找，进行基因注释和功能分析。”
---

# 基因数据库

## 概述

NCBI Gene 是一个整合了不同物种基因信息的综合数据库。它提供命名法、参考序列 (RefSeqs)、染色体图谱、生物途径、遗传变异、表型以及全球基因组资源的交叉引用。

## 何时使用此技能

在处理基因数据时应使用此技能，包括按基因符号或 ID 搜索、检索基因序列和元数据、分析基因功能和通路或执行批量基因查找。

## 快速入门

NCBI 提供两个主要的基因数据访问 API：

1. **电子实用程序**（传统）：适用于所有 Entrez 数据库的全功能 API，具有灵活的查询功能
2. **NCBI 数据集 API**（较新）：通过简化的工作流程针对基因数据检索进行了优化

选择电子实用程序进行复杂查询和跨数据库搜索。选择数据集 API，通过单个请求中的元数据和序列进行简单的基因数据检索。

## 常见工作流程

### 按符号或名称搜索基因

要在生物体中按符号或名称搜索基因：

1. 将 `scripts/query_gene.py` 脚本与 E-utilities ESearch 结合使用
2. 指定基因符号和生物体（例如，“BRCA1 in human”）
3. 脚本返回匹配的基因ID

查询模式示例：
- 基因符号：`insulin[gene name] AND human[organism]`
- 患有疾病的基因：`dystrophin[gene name] AND muscular dystrophy[disease]`
- 染色体位置：`human[organism] AND 17q21[chromosome]`

### 通过ID检索基因信息

要获取已知基因 ID 的详细信息：

1. 将 `scripts/fetch_gene_data.py` 与 Datasets API 结合使用以获取全面的数据
2. 或者，将 `scripts/query_gene.py` 与电子实用程序 EFetch 结合使用以获取特定格式
3. 指定所需的输出格式（JSON、XML 或文本）

数据集 API 返回：
- 基因命名法和别名
- 转录本和蛋白质的参考序列 (RefSeqs)
- 染色体定位和绘图
- 基因本体（GO）注释
- 相关出版物

### 批量基因查找

同时对于多个基因：

1.使用`scripts/batch_gene_lookup.py`进行高效的批处理
2. 提供基因符号或ID列表
3. 为基于符号的查询指定有机体
4. 脚本自动处理速率限制（使用 API 密钥每秒 10 个请求）

此工作流程适用于：
- 验证基因列表
- 检索基因面板的元数据
- 交叉引用基因标识符
- 建立基因注释表

### 按生物背景搜索

寻找与特定生物功能或表型相关的基因：

1. 将电子实用程序与基因本体 (GO) 术语或表型关键字结合使用
2. 通过通路名称或疾病关联查询
3. 按生物体、染色体或其他属性过滤

搜索示例：
- 按 GO 术语：`GO:0006915[biological process]`（细胞凋亡）
- 按表型：`diabetes[phenotype] AND mouse[organism]`
- 按途径：`insulin signaling pathway[pathway]`

### API 访问模式

**速率限制：**
- 没有 API 密钥：电子实用程序每秒 3 个请求，数据集 API 每秒 5 个请求
- 使用 API 密钥：两个 API 每秒 10 个请求

**身份验证：**
在 https://www.ncbi.nlm.nih.gov/account/ 注册免费 NCBI API 密钥以提高速率限制。

**错误处理：**
两个 API 都返回标准 HTTP 状态代码。常见错误包括：
- 400：查询格式错误或参数无效
- 429：超出速率限制
- 404：未找到基因 ID

使用指数退避重试失败的请求。

## 脚本使用

### query_gene.py

使用电子实用程序（ESearch、ESummary、EFetch）查询 NCBI 基因。

```bash
python scripts/query_gene.py --search "BRCA1" --organism "human"
python scripts/query_gene.py --id 672 --format json
python scripts/query_gene.py --search "insulin[gene] AND diabetes[disease]"
```

### fetch_gene_data.py

使用 NCBI 数据集 API 获取全面的基因数据。

<<<代码块_1>>>

###batch_gene_lookup.py

有效处理多个基因查询。

<<<代码块_2>>>

## API 参考

有关详细的 API 文档，包括端点、参数、响应格式和示例，请参阅：

- `references/api_reference.md` - 电子实用程序和数据集 API 的综合 API 文档
- `references/common_workflows.md` - 其他示例和用例模式

当需要特定 API 端点详细信息、参数选项或响应结构信息时，请搜索这些参考。

## 数据格式

NCBI 基因数据可以多种格式检索：

- **JSON**：适合编程处理的结构化数据
- **XML**：具有完整元数据的详细分层格式
- **GenBank**：带注释的序列数据
- **FASTA**：仅序列数据
- **文本**：人类可读的摘要

为现代应用程序选择 JSON，为需要详细元数据的遗留系统选择 XML，为序列分析工作流程选择 FASTA。

## 最佳实践

1. **在通过基因符号搜索时始终指定生物体**以避免歧义
2. **使用基因 ID** 在可用时进行精确查找
3. **处理多个基因时的批量请求**，以最大限度地减少 API 调用
4. **在本地缓存结果**，减少冗余查询
5. **在脚本中包含 API 密钥**以获得更高的速率限制
6. **通过针对暂时性故障的重试逻辑，优雅地处理错误**
7. **在批处理之前验证基因符号**以捕获拼写错误

## 资源

该技能包括：

### 脚本/
- `query_gene.py` - 使用电子实用程序（ESearch、ESummary、EFetch）查询基因
- `fetch_gene_data.py` - 使用 NCBI 数据集 API 获取基因数据
- `batch_gene_lookup.py` - 高效处理多个基因查询

###参考资料/
- `api_reference.md` - 电子实用程序和数据集 API 的详细 API 文档
- `common_workflows.md` - 常见基因查询和用例示例