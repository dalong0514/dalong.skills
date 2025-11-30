<!-- 此文件由机器翻译自 id_mapping_databases.md -->

# UniProt ID 映射数据库

UniProt ID 映射服务支持的数据库的完整列表。调用 ID 映射 API 时使用这些数据库名称。

## 以编程方式检索数据库列表

```python
import requests
response = requests.get("https://rest.uniprot.org/configure/idmapping/fields")
databases = response.json()
```

## UniProt 数据库

### UniProtKB
- `UniProtKB_AC-ID` - UniProt 加入和 ID
- `UniProtKB` - UniProt 知识库
- `UniProtKB-Swiss-Prot` - 已审核（Swiss-Prot）
- `UniProtKB-TrEMBL` - 未经审核 (TrEMBL)
- `UniParc` - UniProt 存档
- `UniRef50` - UniRef 50% 身份集群
- `UniRef90` - UniRef 90% 身份集群
- `UniRef100` - UniRef 100% 身份集群

## 序列数据库

### 核苷酸序列
- `EMBL` - EMBL/GenBank/DDBJ
- `EMBL-CDS` - EMBL 编码序列
- `RefSeq_Nucleotide` - RefSeq 核苷酸序列
- `CCDS` - 共识 CDS

### 蛋白质序列
- `RefSeq_Protein` - RefSeq 蛋白质序列
- `PIR` - 蛋白质信息资源

## 基因数据库

- `GeneID` - Entrez 基因
- `Gene_Name` - 基因名称
- `Gene_Synonym` - 基因同义词
- `Gene_OrderedLocusName` - 有序轨迹名称
- `Gene_ORFName` - ORF 名称

## 基因组数据库

### 一般
- `Ensembl` - 整体
- `EnsemblGenomes` - 整体基因组
- `EnsemblGenomes_PRO` - Ensembl 基因组蛋白质
- `EnsemblGenomes_TRS` - Ensembl 基因组转录本
- `Ensembl_PRO` - 整体蛋白质
- `Ensembl_TRS` - 整体转录

### 特定生物体
- `KEGG` - KEGG 基因
- `PATRIC` - 帕特里克
- `UCSC` - UCSC 基因组浏览器
- `VectorBase` - VectorBase
- `WBParaSite` - WormBase ParaSite

## 结构数据库

- `PDB` - 蛋白质数据库
- `AlphaFoldDB` - AlphaFold 数据库
- `BMRB` - 生物磁共振数据库
- `PDBsum` - PDB 摘要
- `SASBDB` - 小角散射生物数据库
- `SMR` - SWISS-MODEL 存储库

## 蛋白质家族和结构域数据库

- `InterPro` - InterPro
- `Pfam` - Pfam 蛋白家族
- `PROSITE` - PROSITE
- `SMART` - 智能域
- `CDD` - 保守域数据库
- `HAMAP` - HAMAP
- `PANTHER` - 黑豹
- `PRINTS` - 打印
- `ProDom` - ProDom
- `SFLD` - 结构-功能链接数据库
- `SUPFAM` - 超级家庭
- `TIGRFAMs` - TIGRFAM

## 生物体特异性数据库

### 模型生物
- `MGI` - 小鼠基因组信息学
- `RGD` - 大鼠基因组数据库
- `FlyBase` - FlyBase（果蝇）
- `WormBase` - WormBase（线虫）
- `Xenbase` - Xenbase (非洲爪蟾)
- `ZFIN` - 斑马鱼信息网
- `dictyBase` - dictyBase（盘基网柄菌）
- `EcoGene` - EcoGene（大肠杆菌）
- `SGD` - 酵母基因组数据库（酵母）
- `PomBase` - PomBase（粟酒裂殖酵母）
- `TAIR` - 拟南芥信息资源

### 人类特异性
- `HGNC` - HUGO 基因命名委员会
- `CCDS` - 共识编码序列数据库

## 通路数据库

- `Reactome` - 反应组
- `BioCyc` - BioCyc
- `PlantReactome` - 植物反应组
- `SIGNOR` - 签名者
- `SignaLink` - SignaLink

## 酶和新陈代谢

- `EC` - 酶佣金编号
- `BRENDA` - BRENDA 酶数据库
- `SABIO-RK` - SABIO-RK（生化反应）
- `MetaCyc` - MetaCyc

## 疾病和表型数据库

- `OMIM` - 人类在线孟德尔遗传
- `MIM` - MIM（与 OMIM 相同）
- `OrphaNet` - Orphanet（罕见疾病）
- `DisGeNET` - DisGeNET
- `MalaCards` - MalaCards
- `CTD` - 比较毒理学数据库
- `OpenTargets` - 打开目标

## 药物和化学数据库

- `ChEMBL` - ChEMBL
- `DrugBank` - DrugBank
- `DrugCentral` - DrugCentral
- `GuidetoPHARMACOLOGY` - 药理学指南
- `SwissLipids` - SwissLipids

## 基因表达数据库

- `Bgee` - Bgee 基因表达
- `ExpressionAtlas` - 表达式图集
- `Genevisible` - Genevisible
- `CleanEx` - CleanEx
## 蛋白质组数据库

- `PRIDE` - PRIDE 蛋白质组学
- `PeptideAtlas` - 肽图谱
- `ProteomicsDB` - ProteomicsDB
- `CPTAC` - CPTAC
- `jPOST` - jPOST
- `MassIVE` - MassIVE
- `MaxQB` - MaxQB
- `PaxDb` - PaxDb
- `TopDownProteomics` - 自上而下的蛋白质组学

## 蛋白质-蛋白质相互作用

- `STRING` - STRING
- `BioGRID` - BioGRID
- `IntAct` - 完整
- `MINT` - 完好
- `DIP` - 相互作用蛋白质数据库
- `ComplexPortal` - 复杂门户

## 本体论

- `GO` - 基因本体
- `GeneTree` - Ensembl 基因树
- `HOGENOM` - HOGENOM
- `HOVERGEN` - HOVERGEN
- `KO` - KEGG 直系学
- `OMA` - OMA 正系学
- `OrthoDB` - OrthoDB
- `TreeFam` - TreeFam

## 其他专业数据库

### 糖基化
- `GlyConnect` - GlyConnect
- `GlyGen` - GlyGen

### 蛋白质修饰
- `PhosphoSitePlus` - PhosphoSitePlus
- `iPTMnet` - iPTMnet

### 抗体
- `Antibodypedia` - 抗体百科
- `DNASU` - DNASU

### 蛋白质定位
- `COMPARTMENTS` - 隔间
- `NeXtProt` - NeXtProt（人类蛋白质）

### 进化和系统发育
- `eggNOG` - EggNOG
- `GeneTree` - Ensembl 基因树
- `InParanoid` - InParanoid

### 技术资源
- `PRO` - 蛋白质本体论
- `GenomeRNAi` - 基因组RNAi
- `PubMed` - PubMed 文献参考

## 常见的映射场景

### 示例 1：UniProt 到 PDB
<<<代码块_1>>>

### 示例 2：UniProt 的基因名称
<<<代码块_2>>>

### 示例 3：UniProt 到 Ensembl
<<<代码块_3>>>

### 示例 4：RefSeq 到 UniProt
<<<代码块_4>>>

### 示例 5：UniProt 到 GO 术语
<<<代码块_5>>>

## 使用说明

1. **数据库名称区分大小写**：使用列出的确切名称

2. **多对多映射**：一个ID可以映射到多个目标ID

3. **映射失败**：有些ID可能没有映射；检查结果中的 `failedIds` 字段

4. **批量大小限制**：每个作业最多 100,000 个 ID

5. **结果过期**：结果保存7天

6. **双向映射**：大多数数据库支持双向映射

## API 端点

### 获取可用数据库
<<<代码块_6>>>

### 提交测绘作业
```
POST https://rest.uniprot.org/idmapping/run
Content-Type: application/x-www-form-urlencoded

from={from_db}&to={to_db}&ids={comma_separated_ids}
```

### 检查作业状态
```
GET https://rest.uniprot.org/idmapping/status/{jobId}
```

### 获取结果
```
GET https://rest.uniprot.org/idmapping/results/{jobId}
```

## 资源

- ID 映射工具：https://www.uniprot.org/id-mapping
- API 文档：https://www.uniprot.org/help/id_mapping
- 编程访问：https://www.uniprot.org/help/api_idmapping