<!-- 此文件由机器翻译自 api_fields.md -->

# UniProt API 字段参考

用于自定义 UniProt API 查询的可用字段的完整列表。将这些字段与 `fields` 参数结合使用可仅检索您需要的数据。

## 用法

将 fields 参数添加到您的查询中：
```
https://rest.uniprot.org/uniprotkb/search?query=insulin&fields=accession,gene_names,organism_name,length
```

多个字段以逗号分隔。逗号后不能有空格。

## 核心领域

### 身份识别
- `accession` - 主要入藏号（例如，P12345）
- `id` - 条目名称（例如 INSR_HUMAN）
- `uniprotkb_id` - 与 id 相同
- `entryType` - 已审核 (Swiss-Prot) 或未审核 (TrEMBL)

### 蛋白质名称
- `protein_name` - 推荐和替代蛋白质名称
- `gene_names` - 基因名称
- `gene_primary` - 主要基因名称
- `gene_synonym` - 基因同义词
- `gene_oln` - 有序轨迹名称
- `gene_orf` - ORF 名称

### 生物信息
- `organism_name` - 生物体学名
- `organism_id` - NCBI 分类标识符
- `lineage` - 分类谱系
- `virus_hosts` - 病毒宿主生物体（针对病毒蛋白）

### 序列信息
- `sequence` - 氨基酸序列
- `length` - 序列长度
- `mass` - 分子量（道尔顿）
- `fragment` - 条目是否是片段
- `checksum` - 序列 CRC64 校验和

## 注释字段

### 功能和生物学
- `cc_function` - 功能说明
- `cc_catalytic_activity` - 催化活性
- `cc_activity_regulation` - 活动调节
- `cc_pathway` - 代谢途径信息
- `cc_cofactor` - 辅因子信息

### 交互和本地化
- `cc_interaction` - 蛋白质-蛋白质相互作用
- `cc_subunit` - 子单元结构
- `cc_subcellular_location` - 亚细胞位置
- `cc_tissue_specificity` - 组织特异性
- `cc_developmental_stage` - 发育阶段表达

### 疾病和表型
- `cc_disease` - 疾病关联
- `cc_disruption_phenotype` - 破坏表型
- `cc_allergen` - 过敏原信息
- `cc_toxic_dose` - 有毒剂量信息

### 翻译后修饰
- `cc_ptm` - 翻译后修饰
- `cc_mass_spectrometry` - 质谱数据

### 其他评论
- `cc_alternative_products` - 替代产品（亚型）
- `cc_polymorphism` - 多态性信息
- `cc_rna_editing` - RNA 编辑
- `cc_caution` - 注意事项
- `cc_miscellaneous` - 其他信息
- `cc_similarity` - 序列相似性
- `cc_sequence_caution` - 顺序警告
- `cc_web_resource` - 网络资源

## 特征字段 (ft_)

### 分子加工
- `ft_signal` - 信号肽
- `ft_transit` - 转运肽
- `ft_init_met` - 引发剂蛋氨酸
- `ft_propep` - 前肽
- `ft_chain` - 链（成熟蛋白质）
- `ft_peptide` - 肽

### 地区和站点
- `ft_domain` - 域
- `ft_repeat` - 重复
- `ft_ca_bind` - 钙结合
- `ft_zn_fing` - 锌指
- `ft_dna_bind` - DNA 结合
- `ft_np_bind` - 核苷酸结合
- `ft_region` - 感兴趣区域
- `ft_coiled` - 线圈
- `ft_motif` - 短序列基序
- `ft_compbias` - 成分偏差

### 站点和修改
- `ft_act_site` - 活动站点
- `ft_metal` - 金属装订
- `ft_binding` - 结合位点
- `ft_site` - 网站
- `ft_mod_res` - 修饰残基
- `ft_lipid` - 脂化
- `ft_carbohyd` - 糖基化
- `ft_disulfid` - 二硫键
- `ft_crosslnk` - 交叉链接

### 结构特点
- `ft_helix` - 螺旋
- `ft_strand` - Beta 链
- `ft_turn` - 转向
- `ft_transmem` - 跨膜区域
- `ft_intramem` - 膜内区域
- `ft_topo_dom` - 拓扑域

### 变化和冲突
- `ft_variant` - 自然变体
- `ft_var_seq` - 替代序列
- `ft_mutagen` - 诱变
- `ft_unsure` - 不确定残留
- `ft_conflict` - 序列冲突
- `ft_non_cons` - 非连续残基
- `ft_non_ter` - 非末端残基
- `ft_non_std` - 非标准残留

## 基因本体论（GO）

- `go` - 所有 GO 术语
- `go_p` - 生物过程
- `go_c` - 蜂窝组件
- `go_f` - 分子功能
- `go_id` - GO 术语标识符

## 交叉引用 (xref_)

### 序列数据库
- `xref_embl` - EMBL/GenBank/DDBJ
- `xref_refseq` - RefSeq
- `xref_ccds` - CCDS
- `xref_pir` - PIR

### 3D 结构数据库
- `xref_pdb` - 蛋白质数据库
- `xref_pcddb` - PCD 数据库
- `xref_alphafolddb` - AlphaFold 数据库
- `xref_smr` - SWISS-MODEL 存储库

### 蛋白质家族/域数据库
- `xref_interpro` - InterPro
- `xref_pfam` - Pfam
- `xref_prosite` - PROSITE
- `xref_smart` - 智能

### 基因组数据库
- `xref_ensembl` - 整体
- `xref_ensemblgenomes` - 整体基因组
- `xref_geneid` - Entrez 基因
- `xref_kegg` - KEGG

### 生物体特异性数据库
- `xref_mgi` - MGI（鼠标）
- `xref_rgd` - RGD（大鼠）
- `xref_flybase` - FlyBase（飞行）
- `xref_wormbase` - WormBase（蠕虫）
- `xref_xenbase` - Xenbase（青蛙）
- `xref_zfin` - ZFIN（斑马鱼）

### 通路数据库
- `xref_reactome` - 反应组
- `xref_signor` - 签名者
- `xref_signalink` - SignaLink

### 疾病数据库
- `xref_disgenet` - DisGeNET
- `xref_malacards` - MalaCards
- `xref_omim` - OMIM
- `xref_orphanet` - 孤儿

### 药物数据库
- `xref_chembl` - ChEMBL
- `xref_drugbank` - DrugBank
- `xref_guidetopharmacology` - 药理学指南

### 表达数据库
- `xref_bgee` - Bgee
- `xref_expressionetatlas` - 表达式图集
- `xref_genevisible` - Genevisible

## 元数据字段

### 日期
- `date_created` - 条目创建日期
- `date_modified` - 上次修改日期
- `date_sequence_modified` - 上次序列修改日期

### 证据和质量
- `annotation_score` - 注释分数 (1-5)
- `protein_existence` - 蛋白质存在水平
- `reviewed` - 条目是否经过审核 (Swiss-Prot)

###文学
- `lit_pubmed_id` - PubMed 标识符
- `lit_doi` - DOI 标识符

### 蛋白质组学
- `proteome` - 蛋白质组标识符
- `tools` - 用于注释的工具

## 以编程方式检索可用字段

使用配置端点获取所有可用字段：
<<<代码块_1>>>

或者用Python：
<<<代码块_2>>>

## 常见字段组合

### 蛋白质基本信息
<<<代码块_3>>>

### 顺序和结构
<<<代码块_4>>>

### 函数注释
<<<代码块_5>>>

###疾病信息
<<<代码块_6>>>

### 表达模式
```
fields=accession,gene_names,cc_tissue_specificity,cc_developmental_stage,xref_bgee
```

### 完整注释
```
fields=accession,id,protein_name,gene_names,organism_name,sequence,length,cc_*,ft_*,go,xref_pdb
```

## 注释

1. **通配符**：某些字段支持通配符（例如，`cc_*` 适用于所有注释字段，`ft_*` 适用于所有功能）

2. **性能**：请求更少的字段可以提高响应时间并减少带宽

3. **格式依赖性**：根据输出格式（JSON 与 TSV），某些字段的格式可能会有所不同

4. **空值**：响应中可能会省略没有数据的字段 (JSON) 或为空 (TSV)

5. **数组与字符串**：在 JSON 格式中，许多字段返回对象数组而不是简单的字符串

## 资源

- 交互式现场浏览器：https://www.uniprot.org/api-documentation
- API 字段端点：https://rest.uniprot.org/configure/uniprotkb/result-fields
- 返回字段文档：https://www.uniprot.org/help/return_fields