<!-- 此文件由机器翻译自 string_reference.md -->

# STRING 数据库 API 参考

## 概述

STRING（相互作用基因/蛋白质检索搜索工具）是一个综合数据库，包含已知和预测的蛋白质-蛋白质相互作用，整合了来自 40 多个来源的数据。

**数据库统计（v12.0+）：**
- 覆盖范围：5000+基因组
- 蛋白质：约 5930 万
- 互动次数：20+十亿
- 数据类型：物理交互、功能关联、共表达、共现、文本挖掘、数据库

**核心数据资源：** 全球生物数据联盟和 ELIXIR 指定

## API 基本 URL

- **当前版本**：https://string-db.org/api
- **版本特定**：https://version-12-0.string-db.org/api（为了重现性）
- **API 文档**：https://string-db.org/help/api/

## 最佳实践

1. **标识符映射**：始终首先使用 `get_string_ids` 映射标识符，以便更快地进行后续查询
2. **使用 STRING ID**：优先使用 STRING 标识符（例如 `9606.ENSP00000269305`）而不是基因名称
3. **指定物种**：对于具有 >10 个蛋白质的网络，始终指定物种 NCBI 分类单元 ID
4. **速率限制**：API 调用之间等待 1 秒以避免服务器过载
5. **版本化 URL**：使用特定于版本的 URL 进行可重复的研究
6. **POST over GET**：对大型蛋白质列表使用 POST 请求
7. **来电者身份**：包含用于跟踪的 `caller_identity` 参数（例如，您的应用程序名称）

## API 方法

### 1. 标识符映射 (`get_string_ids`)

**目的**：将常见蛋白质名称、基因符号、UniProt ID 和其他标识符映射到 STRING 标识符。

**端点**：`/api/tsv/get_string_ids`

**参数**：
- `identifiers`（必需）：以换行符分隔的蛋白质名称/ID (`%0d`)
- `species`（必需）：NCBI 分类单元 ID
- `limit`：每个标识符的匹配数（默认值：1）
- `echo_query`：在输出中包含查询项（1 或 0）
- `caller_identity`：应用程序标识符

**输出格式**：带有列的 TSV：
- `queryItem`：原始查询
- `queryIndex`：查询位置
- `stringId`：STRING 标识符
- `ncbiTaxonId`：物种分类单元 ID
- `taxonName`：物种名称
- `preferredName`：首选基因名称
- `annotation`：蛋白质描述

**示例**：
```
identifiers=TP53%0dBRCA1&species=9606&limit=1
```

**用例**：
- 将基因符号转换为字符串 ID
- 验证蛋白质标识符
- 查找规范的蛋白质名称

### 2. 网络数据 (`network`)

**目的**：以表格格式检索蛋白质-蛋白质相互作用网络数据。

**端点**：`/api/tsv/network`

**参数**：
- `identifiers`（必需）：由 `%0d` 分隔的蛋白质 ID
- `species`：NCBI 分类单元 ID
- `required_score`：置信度阈值 0-1000（默认值：400）
  - 150：低置信度
  - 400：中等可信度
  - 700：高置信度
  - 900：最高置信度
- `network_type`：`functional`（默认）或`physical`
- `add_nodes`：添加 N 个相互作用的蛋白质 (0-10)
- `caller_identity`：应用程序标识符

**输出格式**：带有列的 TSV：
- `stringId_A`、`stringId_B`：相互作用的蛋白质
- `preferredName_A`、`preferredName_B`：基因名称
- `ncbiTaxonId`：物种
- `score`：综合交互得分 (0-1000)
- `nscore`：邻里分数
- `fscore`：融合得分
- `pscore`：系统发育谱评分
- `ascore`：共表达分数
- `escore`：实验分数
- `dscore`：数据库分数
- `tscore`：文本挖掘分数

**网络类型**：
- **功能性**：所有交互证据类型（推荐用于大多数分析）
- **物理**：仅直接物理结合证据

**示例**：
<<<代码块_1>>>

### 3. 网络图像 (`image/network`)

**目的**：生成 PNG 图像形式的视觉网络表示。

**端点**：`/api/image/network`

**参数**：
- `identifiers`（必需）：由 `%0d` 分隔的蛋白质 ID
- `species`：NCBI 分类单元 ID
- `required_score`：置信度阈值 0-1000
- `network_flavor`：可视化风格
  - `evidence`：将证据类型显示为彩色线
  - `confidence`：将置信度显示为线条粗细
  - `actions`：显示激活/抑制交互
- `add_nodes`：添加 N 个相互作用的蛋白质 (0-10)
- `caller_identity`：应用程序标识符

**输出**：PNG图像（二进制数据）

**图像规格**：
- 格式：PNG
- 大小：根据网络大小自动缩放
- 提供高分辨率选项（添加`?highres=1`）

**示例**：
<<<代码块_2>>>

### 4. 互动合作伙伴 (`interaction_partners`)

**目的**：检索给定蛋白质的所有 STRING 相互作用伙伴。

**端点**：`/api/tsv/interaction_partners`

**参数**：
- `identifiers`（必需）：蛋白质 ID
- `species`：NCBI 分类单元 ID
- `required_score`：置信度阈值 0-1000
- `limit`：最大合作伙伴数量（默认值：10）
- `caller_identity`：应用程序标识符

**输出格式**：与 `network` 方法具有相同列的 TSV

**用例**：
- 寻找枢纽蛋白
- 扩大网络
- 发现新颖的互动

**示例**：
<<<代码块_3>>>

### 5. 功能丰富 (`enrichment`)

**目的**：对多个注释数据库中的一组蛋白质进行功能富集分析。

**端点**：`/api/tsv/enrichment`

**参数**：
- `identifiers`（必需）：蛋白质 ID 列表
- `species`（必需）：NCBI 分类单元 ID
- `caller_identity`：应用程序标识符

**丰富类别**：
- **基因本体**：生物过程、分子功能、细胞成分
- **KEGG 通路**：代谢和信号通路
- **Pfam**：蛋白质结构域
- **InterPro**：蛋白质家族和结构域
- **SMART**：域架构
- **UniProt 关键字**：精选的功能关键字

**输出格式**：带有列的 TSV：
- `category`：注释类别
- `term`：术语 ID
- `description`：术语描述
- `number_of_genes`：输入中包含该术语的基因
- `number_of_genes_in_background`：包含该术语的基因总数
- `ncbiTaxonId`：物种
- `inputGenes`：逗号分隔的基因列表
- `preferredNames`：逗号分隔的基因名称
- `p_value`：富集 p 值（未校正）
- `fdr`：错误发现率（校正后的 p 值）

**统计方法**：采用 Benjamini-Hochberg FDR 校正的 Fisher 精确检验

**示例**：
<<<代码块_4>>>

### 6. PPI 丰富 (`ppi_enrichment`)

**目的**：测试网络的交互是否比偶然预期的要多得多。

**端点**：`/api/json/ppi_enrichment`

**参数**：
- `identifiers`（必需）：蛋白质 ID 列表
- `species`：NCBI 分类单元 ID
- `required_score`：置信度阈值
- `caller_identity`：应用程序标识符

**输出格式**：带有字段的 JSON：
- `number_of_nodes`：网络中的蛋白质
- `number_of_edges`：观察到的交互
- `expected_number_of_edges`：预期交互（随机）
- `p_value`：统计意义

**解释**：
- p 值 < 0.05：网络显着丰富
- 低 p 值表明蛋白质形成功能模块

**示例**：
<<<代码块_5>>>

### 7. 同源分数 (`homology`)

**目的**：检索蛋白质相似性/同源性得分。

**端点**：`/api/tsv/homology`

**参数**：
- `identifiers`（必需）：蛋白质 ID
- `species`：NCBI 分类单元 ID
- `caller_identity`：应用程序标识符

**输出格式**：TSV 与蛋白质之间的同源性得分

**用例**：
- 识别蛋白质家族
- 旁系同源分析
- 跨物种比较

### 8.版本信息(`version`)

**用途**：返回当前 STRING 数据库版本。

**端点**：`/api/tsv/version`

**输出**：版本字符串（例如“12.0”）

## 常见物种 NCBI 分类单元 ID

|有机体 |通用名称|分类单元 ID |
|----------|-------------|----------|
|智人 |人类 | 9606 | 9606
|小家鼠 |鼠标| 10090 | 10090
|褐家鼠 |老鼠 | 10116 |
|果蝇 |果蝇| 7227 | 7227
|秀丽隐杆线虫 |线虫 | 6239 | 6239
|酿酒酵母|酵母| 4932 |
|拟南芥|塔勒水芹 | 3702 | 3702
|大肠杆菌 K-12 |大肠杆菌 | 511145 | 511145
|斑马鱼 |斑马鱼 | 7955 | 7955
|鸡内金 |鸡 | 9031 | 9031

完整列表：https://string-db.org/cgi/input?input_page_active_form=organisms

## STRING 标识符格式

STRING 使用带有分类单元前缀的 Ensembl 蛋白质 ID：
- 格式：`{taxonId}.{ensemblProteinId}`
- 示例：`9606.ENSP00000269305`（人类 TP53）

**ID 组件**：
- **分类单元 ID**：NCBI 分类标识符
- **蛋白质 ID**：通常为 Ensembl 蛋白质 ID (ENSP...)

## 交互置信度分数

STRING 提供基于多个证据渠道的综合置信度分数 (0-1000)：

### 证据渠道
1. **邻域（nscore）**：基因融合和保守的基因组邻域
2. **融合（fscore）**：跨物种的基因融合事件
3. **系统发育概况（pscore）**：跨物种共现
4. **共表达（ascore）**：RNA表达相关性
5. **实验（escore）**：生化/遗传实验
6. **数据库 (dscore)**：精心策划的途径/复杂数据库
7. **文本挖掘（tscore）**：文献共现

### 推荐阈值

- **150**：低置信度（探索性分析）
- **400**：中等置信度（标准分析）
- **700**：高置信度（保守分析）
- **900**：最高置信度（非常严格）

## 输出格式

### 可用格式

1. **TSV**：制表符分隔值（默认，最适合数据处理）
2. **JSON**：JavaScript 对象表示法（结构化数据）
3. **XML**：可扩展标记语言
4. **PSI-MI**：蛋白质组学标准倡议格式
5. **PSI-MITAB**：制表符分隔的 PSI-MI 格式
6. **PNG**：图像格式（用于网络可视化）
7. **SVG**：可缩放矢量图形（用于网络可视化）

### 格式选择

将 URL 中的 `/tsv/` 替换为所需格式：
- `/json/network` - JSON 格式
- `/xml/network` - XML 格式
- `/image/network` - PNG 图像

## 错误处理

### HTTP 状态代码

- **200 OK**：请求成功
- **400 错误请求**：无效的参数或语法
- **404 Not Found**：未找到蛋白质/物种
- **500 内部服务器错误**：服务器错误

### 常见错误

1. **“未发现蛋白质”**：标识符无效或物种不匹配
2. **“所需物种”**：大型网络缺少物种参数
3. **空结果**：没有高于分数阈值的交互
4. **超时**：网络太大，减少蛋白质计数

## 高级功能

### 批量网络上传

对于完整的蛋白质组分析：
1. 导航至https://string-db.org/
2. 选择“上传蛋白质组”选项
3.上传FASTA文件
4. STRING生成完整的交互网络并预测功能

### 值/排名丰富 API

对于差异表达/蛋白质组数据：

1. **获取API密钥**：
<<<代码块_6>>>

2. **提交数据**：制表符分隔的蛋白质 ID 和值对

3. **检查状态**：
```
/api/json/valuesranks_enrichment_status?job_id={id}
```

4. **检索结果**：访问丰富的表格和图形

**要求**：
- 完整的蛋白质组（无过滤）
- 每种蛋白质的数值
- 正确的物种标识符

### 网络定制

**网络大小控制**：
- `add_nodes=N`：添加 N 个最相关的蛋白质
- `limit`：控制伙伴检索

**置信过滤**：
- 根据分析目标调整`required_score`
- 分数越高=假阳性越少，假阴性越多

**网络类型选择**：
- `functional`：所有证据（建议用于路径分析）
- `physical`：仅直接结合（推荐用于结构研究）

## 与其他工具集成

### Python 库

**请求**（推荐）：
```python
import requests
url = "https://string-db.org/api/tsv/network"
params = {"identifiers": "TP53", "species": 9606}
response = requests.get(url, params=params)
```

**urllib**（标准库）：
```python
import urllib.request
url = "https://string-db.org/api/tsv/network?identifiers=TP53&species=9606"
response = urllib.request.urlopen(url)
```

### R 集成

**STRINGdb Bioconductor 包**：
```R
library(STRINGdb)
string_db <- STRINGdb$new(version="12", species=9606)
```

### 细胞景观

STRING 网络可以导入 Cytoscape 进行可视化和分析：
1.使用stringApp插件
2.导入TSV网络数据
3.应用布局和样式

## 数据许可

STRING 数据可根据 **Creative Commons BY 4.0** 许可证免费获取：
- ✓ 免费用于学术和商业目的
- ✓ 需要出处
- ✓ 允许修改
- ✓ 允许再分配

**引用**：Szklarczyk 等人。 （最新出版）

## 速率限制和使用

- **速率限制**：没有严格限制，但避免快速请求
- **建议**：通话之间等待 1 秒
- **大型数据集**：使用来自 https://string-db.org/cgi/download 的批量下载
- **蛋白质组规模**：使用网络上传功能而不是 API

## 相关资源

- **STRING 网站**：https://string-db.org
- **下载页面**：https://string-db.org/cgi/download
- **帮助中心**：https://string-db.org/help/
- **API 文档**：https://string-db.org/help/api/
- **出版物**：https://string-db.org/cgi/about

## 故障排除

**没有返回结果**：
- 验证物种参数与标识符匹配
- 检查标识符格式
- 较低的置信阈值
- 首先使用标识符映射

**超时错误**：
- 减少输入蛋白质的数量
- 将大型查询分成批次
- 使用批量下载进行蛋白质组规模分析

**版本不一致**：
- 使用特定于版本的 URL
- 使用 `/version` 端点检查 STRING 版本
- 如果使用旧 ID，则更新标识符