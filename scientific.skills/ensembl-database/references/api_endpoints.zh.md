<!-- 此文件由机器翻译自 api_endpoints.md -->

# Ensembl REST API 端点参考

Ensembl REST API（2025 年 9 月第 115 版）中提供的所有 17 个 API 端点类别的综合文档。

**基本网址：**
- 当前程序集：`https://rest.ensembl.org`
- GRCh37/hg19（人类）：`https://grch37.rest.ensembl.org`

**速率限制：**
- 匿名：15 个请求/秒
- 已注册：55,000 个请求/小时

## 1. 存档

检索有关已退役 Ensembl 标识符的历史信息。

**获取 /archive/id/:id**
- 检索已退役标识符的存档条目
- 示例：`/archive/id/ENSG00000157764`（退役基因 ID）

## 2. 比较基因组学

访问跨物种的基因树、基因组比对和同源性数据。

**获取/对齐/区域/：物种/：区域**
- 获取某个区域的基因组比对
- 示例：`/alignment/region/human/2:106040000-106040050:1?species_set_group=mammals`

**获取/genetree/id/:id**
- 检索基因家族的基因树
- 示例：`/genetree/id/ENSGT00390000003602`

**获取/genetree/member/id/:id**
- 通过成员基因ID获取基因树
- 示例：`/genetree/member/id/ENSG00000139618`

**GET /同源/id/:id**
- 查找基因的直向同源物和旁系同源物
- 参数：`target_species`、`type`（直向同源物、旁系同源物、全部）
- 示例：`/homology/id/ENSG00000139618?target_species=mouse`

**获取/同源/符号/：物种/：符号**
- 通过基因符号查找同源物
- 示例：`/homology/symbol/human/BRCA2?target_species=mouse`

## 3.交叉引用

将外部数据库标识符链接到 Ensebl 对象。

**获取 /xrefs/id/:id**
- 获取 Ensembl ID 的外部参考
- 示例：`/xrefs/id/ENSG00000139618`

**获取 /xrefs/symbol/:species/:symbol**
- 通过基因符号获取交叉引用
- 示例：`/xrefs/symbol/human/BRCA2`

**获取/外部参照/名称/：物种/：名称**
- 按外部名称搜索对象
- 示例：`/xrefs/name/human/NP_000050`

## 4. 信息

查询有关物种、组件、生物型和数据库版本的元数据。

**获取/信息/物种**
- 列出所有可用物种
- 返回物种名称、组件、分类 ID

**获取/信息/组件/：物种**
- 获取物种的组装信息
- 示例：`/info/assembly/human`（返回 GRCh38.p14）

**获取/信息/组件/：物种/：区域**
- 获取有关染色体区域的详细信息
- 示例：`/info/assembly/human/X`

**获取/信息/生物型/：物种**
- 列出所有可用的生物型（基因类型）
- 示例：`/info/biotypes/human`

**获取/信息/分析/：物种**
- 列出可用的分析类型
- 示例：`/info/analysis/human`

**获取/信息/数据**
- 获取有关当前 Ensebl 版本的一般信息

## 5.连锁不平衡（LD）

计算变异之间的连锁不平衡。

**GET /ld/:species/:id/:population_name**
- 计算变体的 LD
- 示例：`/ld/human/rs1042522/1000GENOMES:phase_3:KHV`

**GET /ld/pairwise/:species/:id1/:id2**
- 计算两个变体之间的LD
- 示例：`/ld/pairwise/human/rs1042522/rs11540652`

## 6. 查找

识别标识符的物种和数据库信息。

**获取/查找/id/:id**
- 通过 Ensembl ID 查找对象
- 参数：`expand`（包括子对象）
- 示例：`/lookup/id/ENSG00000139618?expand=1`

**发布/查找/id**
- 批量查找多个ID
- 提交 ID 的 JSON 数组
- 示例：`{"ids": ["ENSG00000139618", "ENSG00000157764"]}`

**获取/查找/符号/：物种/：符号**
- 按符号查找基因
- 参数：`expand`（包括成绩单）
- 示例：`/lookup/symbol/human/BRCA2?expand=1`

## 7. 映射

转换组件、cDNA、CDS 和蛋白质位置之间的坐标。

**获取/map/cdna/:id/:区域**
- 将 cDNA 坐标映射到基因组
- 示例：`/map/cdna/ENST00000288602/100..300`

**获取 /map/cds/:id/:region**
- 将 CDS 坐标映射到基因组
- 示例：`/map/cds/ENST00000288602/1..300`

**获取/地图/翻译/：id/：区域**
- 将蛋白质坐标映射到基因组
- 示例：`/map/translation/ENSP00000288602/1..100`

**GET /map/:species/:asm_one/:region/:asm_two**
- 装配体之间的地图坐标
- 示例：`/map/human/GRCh37/7:140453136..140453136/GRCh38`

**POST /map/:species/:asm_one/:asm_two**
- 批量装配映射
- 提交区域的 JSON 数组

## 8.本体论和分类法

搜索生物本体和分类学分类。

**获取/本体/id/:id**
- 获取本体术语信息
- 示例：`/ontology/id/GO:0005515`

**GET /本体/名称/:名称**
- 按术语名称搜索本体
- 示例：`/ontology/name/protein%20binding`

**获取/分类/分类/：id**
- 获取分类学分类
- 示例：`/taxonomy/classification/9606`（人类）

**获取/分类/id/:id**
- 通过ID获取分类信息
- 示例：`/taxonomy/id/9606`

## 9. 重叠

查找与某个区域重叠的基因组特征。

**获取/重叠/id/：id**
- 获取与基因/转录本重叠的特征
- 参数：`feature`（基因、转录本、cd、外显子、重复序列等）
- 示例：`/overlap/id/ENSG00000139618?feature=transcript`
**获取/重叠/区域/：物种/：区域**
- 获取基因组区域的所有特征
- 参数：`feature`（基因、转录本、变异、调控等）
- 示例：`/overlap/region/human/7:140424943..140624564?feature=gene`

**获取/重叠/翻译/：id**
- 获取蛋白质特征
- 示例：`/overlap/translation/ENSP00000288602`

## 10.表型注释

检索疾病和性状关联。

**获取/表型/加入/：物种/：加入**
- 通过本体加入获取表型
- 示例：`/phenotype/accession/human/EFO:0003767`

**获取/表型/基因/：物种/：基因**
- 获取基因的表型关联
- 示例：`/phenotype/gene/human/ENSG00000139618`

**获取/表型/区域/：物种/：区域**
- 获取基因组区域的表型
- 示例：`/phenotype/region/human/7:140424943-140624564`

**获取/表型/术语/：物种/：术语**
- 按术语搜索表型
- 示例：`/phenotype/term/human/cancer`

## 11. 监管

访问监管特征和结合基序数据。

**GET /监管/物种/：物种/微阵列/：微阵列/：探针**
- 获取微阵列探针信息
- 示例：`/regulatory/species/human/microarray/HumanWG_6_V2/ILMN_1773626`

**获取/物种/：物种/绑定矩阵/：绑定矩阵_id**
- 获取转录因子结合矩阵
- 示例：`/species/human/binding_matrix/ENSPFM0001`

## 12. 顺序

检索基因组、转录本和蛋白质序列。

**获取/序列/id/:id**
- 通过ID获取序列
- 参数：`type`（基因组、cds、cdna、蛋白质）、`format`（json、fasta、文本）
- 示例：`/sequence/id/ENSG00000139618?type=genomic`

**POST /序列/id**
- 批量序列检索
- 示例：`{"ids": ["ENSG00000139618", "ENSG00000157764"]}`

**获取/序列/区域/：物种/：区域**
- 获取区域的基因组序列
- 参数：`coord_system`、`format`
- 示例：`/sequence/region/human/7:140424943..140624564?format=fasta`

**POST /序列/区域/：物种**
- 批量区域序列检索

## 13. 转录单倍型

从阶段性基因型计算转录单倍型。

**GET /transcript_haplotypes/:species/:id**
- 获取转录单倍型
- 示例：`/transcript_haplotypes/human/ENST00000288602`

## 14.变异效应预测器（VEP）

预测变异的功能后果。

**获取 /vep/:species/hgvs/:hgvs_notation**
- 使用 HGVS 表示法预测变异效应
- 参数：众多 VEP 选项
- 示例：`/vep/human/hgvs/ENST00000288602:c.803C>T`

**POST /vep/:species/hgvs**
- 使用 HGVS 进行批量 VEP 分析
- 示例：`{"hgvs_notations": ["ENST00000288602:c.803C>T"]}`

**获取/vep/:物种/id/:id**
- 预测变体 ID 的影响
- 示例：`/vep/human/id/rs699`

**POST /vep/:物种/id**
- 按变体 ID 批量 VEP

**获取/vep/：物种/区域/：区域/：等位基因**
- 预测区域和等位基因的影响
- 示例：`/vep/human/region/7:140453136:C/T`

**POST /vep/:物种/地区**
- 按地区批量VEP

## 15. 变化

查询遗传变异数据和相关出版物。

**获取/变异/：物种/：id**
- 通过ID获取变体信息
- 参数：`pops`（包括总体频率）、`genotypes`
- 示例：`/variation/human/rs699?pops=1`

**POST /变异/：物种**
- 批量变体查询
- 示例：`{"ids": ["rs699", "rs6025"]}`

**GET /variation/:species/pmcid/:pmcid**
- 从 PubMed Central 文章获取变体
- 示例：`/variation/human/pmcid/PMC5002951`

**GET /variation/:species/pmid/:pmid**
- 从 PubMed 文章中获取变体
- 示例：`/variation/human/pmid/26318936`

## 16. 变体 GA4GH

使用 GA4GH 标准访问基因组变异数据。

**发布/ga4gh/信标**
- 查询信标是否存在变体

**获取/ga4gh/features/:id**
- 通过 ID 获取 GA4GH 格式的特征

**发布/ga4gh/功能/搜索**
- 使用GA4GH协议搜索功能

**POST /ga4gh/variants/search**
- 使用 GA4GH 协议搜索变体

## 响应格式

大多数端点支持多种响应格式：
- **JSON**（默认）：`Content-Type: application/json`
- **FASTA**：用于序列数据
- **XML**：某些端点支持 XML
- **文本**：纯文本输出

使用以下方式指定格式：
1. `Content-Type` 标头
2. URL参数：`content-type=text/x-fasta`
3. 文件扩展名：`/sequence/id/ENSG00000139618.fasta`

## 常用参数

许多端点共享这些参数：

- **展开**：包括子对象（转录本、蛋白质）
- **格式**：输出格式（json、xml、fasta）
- **db_type**：数据库类型（核心、其他功能、变体）
- **object_type**：要返回的对象类型
- **物种**：物种名称（可以是常见的或科学的）

## 错误代码

- **200**：成功
- **400**：错误请求（无效参数）
- **404**：未找到（ID 不存在）
- **429**：超出速率限制
- **500**：内部服务器错误

## 最佳实践

1. **使用批量端点**进行多个查询（更高效）
2. **缓存响应**以最大限度地减少 API 调用
3. **检查响应中的速率限制标头**
4. **通过遵守 `Retry-After` 标头处理 429 错误**
5. **对序列数据使用适当的内容类型**
6. **查询旧基因组版本时指定程序集**
7. **当您需要完整的对象详细信息时启用扩展参数**