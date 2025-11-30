<!-- 此文件由机器翻译自 module_reference.md -->

# gget 模块参考

所有 gget 模块的综合参数参考。

## 参考和基因信息模块

### gget 参考
检索 Ensembl 参考基因组 FTP 和元数据。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `species` | STR | Genus_species 格式的物种或快捷方式（“人类”、“小鼠”）|必填|
| `-w/--which` | STR |要返回的文件类型：gtf、cdna、dna、cds、cdrna、pep |全部 |
| `-r/--release` |整数 |合奏发行号 |最新 |
| `-od/--out_dir` | STR |输出目录路径 |无 |
| `-o/--out` | STR |结果的 JSON 文件路径 |无 |
| `-l/--list_species` |旗帜|列出可用的脊椎动物物种 |假 |
| `-liv/--list_iv_species` |旗帜|列出可用的无脊椎动物物种 |假 |
| `-ftp` |旗帜|仅返回 FTP 链接 |假 |
| `-d/--download` |旗帜|下载文件（需要curl） |假 |
| `-q/--quiet` |旗帜|抑制进度信息 |假 |

**返回：** 包含 FTP 链接、Ensembl 发行号、发行日期、文件大小的 JSON

---

### gget 搜索
在 Ensembl 中按名称或描述搜索基因。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `searchwords` |字符串/列表 |搜索词（不区分大小写）|必填|
| `-s/--species` | STR |目标物种或核心数据库名称|必填|
| `-r/--release` |整数 |合奏发行号 |最新 |
| `-t/--id_type` | STR |返回“基因”或“转录本”| '基因' |
| `-ao/--andor` | STR | “或”（任何术语）或“和”（所有术语）| '或' |
| `-l/--limit` |整数 |返回的最大结果 |无 |
| `-o/--out` | STR |输出文件路径（CSV/JSON）|无 |

**返回：** ensembl_id、gene_name、ensembl_description、ext_ref_description、生物型、URL

---

### 获取信息
从 Ensembl、UniProt 和 NCBI 获取全面的基因/转录本元数据。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `ens_ids` |字符串/列表 | Ensembl ID（也支持 WormBase、Flybase）|必填 |
| `-o/--out` | STR |输出文件路径（CSV/JSON）|无 |
| `-n/--ncbi` |布尔 |禁用 NCBI 数据检索 |假 |
| `-u/--uniprot` |布尔 |禁用 UniProt 数据检索 |假|
| `-pdb` |布尔 |包括 PDB 标识符 |假 |
| `-csv` |旗帜|返回 CSV 格式 (CLI) |假 |
| `-q/--quiet` |旗帜|抑制进度显示 |假 |

**特定于 Python：**
- `save=True`：将输出保存到当前目录
- `wrap_text=True`：使用换行文本格式化数据框

**注意：** 同时处理超过 1000 个 ID 可能会导致服务器错误。

**返回：** UniProt ID、NCBI 基因 ID、基因名称、同义词、蛋白质名称、描述、生物型、规范转录本

---

### gget seq
以 FASTA 格式检索核苷酸或氨基酸序列。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `ens_ids` |字符串/列表 |整体标识符|必填|
| `-o/--out` | STR |输出文件路径 |标准输出|
| `-t/--translate` |旗帜|获取氨基酸序列|假 |
| `-iso/--isoforms` |旗帜|返回所有转录变体 |假 |
| `-q/--quiet` |旗帜|抑制进度信息 |假 |

**数据来源：** Ensembl（核苷酸）、UniProt（氨基酸）

**返回：** FASTA 格式序列

---

## 序列分析和比对模块

### gget爆炸
针对标准数据库的 BLAST 序列。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `sequence` | STR | FASTA/.txt 的序列或路径 |必填|
| `-p/--program` | STR |爆炸n、爆炸p、爆炸x、tblastn、tblastx |自动检测 |
| `-db/--database` | STR | nt、refseq_rna、pdbnt、nr、swissprot、pdbaa、refseq_蛋白质 | nt 或 nr |
| `-l/--limit` |整数 |返回的最大点击数 | 50 | 50
| `-e/--expect` |浮动| E值截止| 10.0 |
| `-lcf/--low_comp_filt` |旗帜|启用低复杂度过滤 |假 |
| `-mbo/--megablast_off` |旗帜|禁用 MegaBLAST（仅限blastn）|假 |
| `-o/--out` | STR |输出文件路径 |无 |
| `-q/--quiet` |旗帜|抑制进步|假|
**返回：** 描述、科学名称、通用名称、Taxid、最高分数、总分、查询覆盖率

---

### gget blat
使用 UCSC BLAT 查找基因组位置。

**参数：**
|参数|类型 |描述 |默认|
|------------|------|-------------|---------|
| `sequence` | STR | FASTA/.txt 的序列或路径 |必填 |
| `-st/--seqtype` | STR | 'DNA'、'蛋白质'、'翻译%20RNA'、'翻译%20DNA' |自动检测 |
| `-a/--assembly` | STR |靶材组装（hg38、mm39、taeGut2等）| '人类'/hg38 |
| `-o/--out` | STR |输出文件路径 |无 |
| `-csv` |旗帜|返回 CSV 格式 (CLI) |假 |
| `-q/--quiet` |旗帜|抑制进步|假|

**返回：**基因组、查询大小、比对开始/结束、匹配、不匹配、比对百分比

---

### 肌肉
使用 Muscle5 比对多个序列。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `fasta` |字符串/列表 |序列或 FASTA 文件路径 |必填|
| `-o/--out` | STR |输出文件路径 |标准输出|
| `-s5/--super5` |旗帜|使用Super5算法（更快，大数据集）|假|
| `-q/--quiet` |旗帜|抑制进步|假 |

**返回：** ClustalW 格式对齐或对齐的 FASTA (.afa)

---

### 获取钻石
快速本地蛋白质/翻译 DNA 比对。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `query` |字符串/列表 |查询序列或FASTA文件|必填 |
| `--reference` |字符串/列表 |参考序列或 FASTA 文件 |必填 |
| `--sensitivity` | STR |快速、中灵敏、灵敏、较灵敏、非常灵敏、超灵敏 |非常敏感|
| `--threads` |整数 | CPU 线程 | 1 |
| `--diamond_binary` | STR | DIAMOND 安装路径 |自动检测 |
| `--diamond_db` | STR |保存数据库以供重复使用 |无 |
| `--translated` |旗帜|实现核苷酸与氨基酸的比对 |假|
| `-o/--out` | STR |输出文件路径 |无 |
| `-csv` |旗帜| CSV 格式 (CLI) |假 |
| `-q/--quiet` |旗帜|抑制进步 |假 |

**返回：** 同一性%、序列长度、匹配位置、空位、E 值、位分数

---

## 结构和蛋白质分析模块

### gget pdb
查询RCSB蛋白质数据库。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `pdb_id` | STR | PDB 标识符（例如“7S7U”）|必填 |
| `-r/--resource` | STR | pdb、条目、pubmed、程序集、实体类型 | 'pdb'|
| `-i/--identifier` | STR |程序集、实体或链 ID |无 |
| `-o/--out` | STR |输出文件路径 |标准输出|

**返回：** PDB 格式（结构）或 JSON（元数据）

---

### gget alphafold
使用 AlphaFold2 预测 3D 蛋白质结构。

**设置：** 需要 OpenMM 和 `gget setup alphafold`（~4GB 下载）

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `sequence` |字符串/列表 |氨基酸序列或 FASTA 文件 |必填 |
| `-mr/--multimer_recycles` |整数 |多聚体的回收迭代| 3 |
| `-o/--out` | STR |输出文件夹路径|带时间戳 |
| `-mfm/--multimer_for_monomer` |旗帜|将多聚体模型应用于单体 |假 |
| `-r/--relax` |旗帜|顶级模特的琥珀放松 |假 |
| `-q/--quiet` |旗帜|抑制进步 |假 |

**仅限 Python：**
- `plot` (bool)：生成 3D 可视化（默认值：True）
- `show_sidechains`（布尔）：包括侧链（默认值：True）

**注意：** 多个序列自动触发多聚体建模

**返回：** PDB结构文件，JSON对齐错误数据，可选3D图

---

### gget 榆树
预测真核线性基序。

**设置：** 需要 `gget setup elm`

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `sequence` | STR |氨基酸序列或 UniProt Acc |必填|
| `-s/--sensitivity` | STR |钻石对准灵敏度 |非常敏感|
| `-t/--threads` |整数 |线程数 | 1 |
| `-bin/--diamond_binary` | STR | DIAMOND 二进制文件的路径 |自动检测 |
| `-o/--out` | STR |输出目录路径|无 |
| `-u/--uniprot` |旗帜|输入是 UniProt Acc |假 |
| `-e/--expand` |旗帜|包括蛋白质名称、生物体、参考文献 |假|
| `-csv` |旗帜| CSV 格式 (CLI) |假|
| `-q/--quiet` |旗帜|抑制进步|假|

**返回：** 两个输出：
1. **ortholog_df**：来自直系同源蛋白的基序
2. **regex_df**：输入序列中匹配的主题

---

## 表达和疾病数据模块

### gget archs4
查询 ARCHS4 的基因相关性或组织表达。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `gene` | STR |基因符号或 Ensembl ID |必填|
| `-w/--which` | STR | “相关性”或“组织”| '相关性' |
| `-s/--species` | STR | “人类”或“小鼠”（仅组织）| '人类' |
| `-o/--out` | STR |输出文件路径 |无 |
| `-e/--ensembl` |旗帜|输入是 Ensembl ID |假|
| `-csv` |旗帜| CSV 格式 (CLI) |假|
| `-q/--quiet` |旗帜|抑制进步|假|

**退货：**
- **相关性**：基因符号、皮尔逊相关系数（前 100 名）
- **组织**：组织 ID、最小/Q1/中值/Q3/最大表达

---

### gget cellxgene
查询 CZ CELLxGENE 发现单细胞数据的普查。

**设置：** 需要 `gget setup cellxgene`

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `--gene` (-g) |列表 |基因名称或 Ensembl ID（区分大小写！）|必填|
| `--tissue` |列表 |组织类型 |无 |
| `--cell_type` |列表 |细胞类型 |无 |
| `--species` (-s) | STR | 'homo_sapiens' 或 'mus_musculus' | '智人' |
| `--census_version`（-cv）| STR | “稳定”、“最新”或过时版本 | “稳定”|
| `-o/--out` | STR |输出文件路径（CLI 必需）|必填|
| `--ensembl` (-e) |旗帜|使用 Ensebl ID |假|
| `--meta_only`（-mo）|旗帜|仅返回元数据 |假|
| `-q/--quiet` |旗帜|抑制进步|假|

**附加过滤器：**疾病、发育阶段、性别、检测、数据集 ID、供体 ID、种族、悬浮类型

**重要提示：** 基因符号区分大小写（“PAX7”表示人类，“Pax7”表示小鼠）

**返回：** 带有计数矩阵和元数据的 AnnData 对象

---

### gget 丰富器
使用 Enrichr/modEnrichr 执行富集分析。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `genes` |列表 |基因符号或 Ensembl ID |必填|
| `-db/--database` | STR |参考数据库或快捷方式 |必填|
| `-s/--species` | STR |人类、小鼠、苍蝇、酵母、蠕虫、鱼 | '人类' |
| `-bkg_l/--background_list` |列表 |背景基因|无 |
| `-o/--out` | STR |输出文件路径 |无 |
| `-ko/--kegg_out` | STR | KEGG 通路图像目录 |无 |

**仅限 Python：**
- `plot` (bool): 生成图形结果

**数据库快捷方式：**
-“途径”→KEGG_2021_Human
-“转录”→ ChEA_2016
-“本体论”→ GO_Biological_Process_2021
-“疾病_药物”→ GWAS_Catalog_2019
-“细胞类型”→ PanglaoDB_Augmented_2021

**返回：** 路径/功能与调整后的 p 值、重叠基因计数的关联

---

### gget bgee
从 Bgee 检索直系同源和表达。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `ens_id` |字符串/列表 | Ensembl 或 NCBI 基因 ID |必填|
| `-t/--type` | STR | '直向同源物' 或 '表达' | '直系同源物' |
| `-o/--out` | STR |输出文件路径 |无 |
| `-csv` |旗帜| CSV 格式 (CLI) |假|
| `-q/--quiet` |旗帜|抑制进步|假|

**注意：** `type='expression'` 时支持多个 ID

**退货：**
- **直向同源物**：跨物种的基因，包含 ID、名称、分类信息
- **表达**：解剖实体、置信度分数、表达状态

---

### gget opentargets
从 OpenTargets 检索疾病/药物关联。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `ens_id` | STR |整体基因 ID |必填|
| `-r/--resource` | STR |疾病、药物、易处理性、药物遗传学、表达、depmap、相互作用 | '疾病' |
| `-l/--limit` |整数 |最大结果 |无 |
| `-o/--out` | STR |输出文件路径 |无 |
| `-csv` |旗帜| CSV 格式 (CLI) |假|
| `-q/--quiet` |旗帜|抑制进步|假|

**特定于资源的过滤器：**
- 药物：`--filter_disease`
- 药物遗传学：`--filter_drug`
- 表达式/depmap：`--filter_tissue`、`--filter_anat_sys`、`--filter_organ`
- 交互：`--filter_protein_a`、`--filter_protein_b`、`--filter_gene_b`

**返回：** 疾病/药物关联、易处理性、药物遗传学、表达、DepMap、相互作用

---

### gget cbio
绘制 cBioPortal 的癌症基因组学热图。

**子命令：** 搜索、绘图

**搜索参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `keywords` |列表 |搜索词|必填|

**绘图参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `-s/--study_ids` |列表 | cBioPortal 研究 ID |必填|
| `-g/--genes` |列表 |基因名称或 Ensembl ID |必填|
| `-st/--stratification` | STR |组织、cancer_type、cancer_type_detailed、study_id、样本 |无 |
| `-vt/--variation_type` | STR |突变_发生、cna_nonbinary、sv_发生、cna_发生、后果 |无 |
| `-f/--filter` | STR |按列值过滤（例如“study_id:msk_impact_2017”）|无 |
| `-dd/--data_dir` | STR |缓存目录| ./gget_cbio_cache |
| `-fd/--figure_dir` | STR |输出目录 | ./gget_cbio_figures |
| `-t/--title` | STR |自定义图形标题 |无 |
| `-dpi` |整数 |分辨率| 100 | 100
| `-q/--quiet` |旗帜|抑制进步|假|
| `-nc/--no_confirm` |旗帜|跳过下载确认 |假|
| `-sh/--show` |旗帜|在窗口中显示绘图 |假|

**返回：** PNG 热图

---

### gget 宇宙
在 COSMIC 数据库中搜索癌症突变。

**重要提示：** 商业用途的许可费用。需要 COSMIC 帐户。

**查询参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `searchterm` | STR |基因名称、Ensembl ID、突变、样本 ID |必填|
| `-ctp/--cosmic_tsv_path` | STR | COSMIC TSV 文件的路径 |必填|
| `-l/--limit` |整数 |最大结果 | 100 | 100
| `-csv` |旗帜| CSV 格式 (CLI) |假|

**下载参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `-d/--download_cosmic` |旗帜|激活下载模式 |假|
| `-gm/--gget_mutate` |旗帜|为 gget mutate 创建版本 |假|
| `-cp/--cosmic_project` | STR |癌症、普查、细胞系、耐药性、基因组筛选、目标筛选 |无 |
| `-cv/--cosmic_version` | STR |宇宙版 |最新 |
| `-gv/--grch_version` |整数 |人类参考基因组（37 或 38）|无 |
| `--email` | STR | COSMIC 账户邮箱 |必填|
| `--password` | STR | COSMIC账户密码 |必填|

**注意：**首次用户必须下载数据库

**返回：** 来自 COSMIC 的突变数据

---

## 附加工具

### gget 变异
生成突变的核苷酸序列。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `sequences` |字符串/列表 | FASTA 文件或序列 |必填|
| `-m/--mutations` | str/df | CSV/TSV 文件或 DataFrame |必填|
| `-mc/--mut_column` | STR |突变列名称 | '突变' |
| `-sic/--seq_id_column` | STR |序列 ID 列 | 'seq_ID' | 'seq_ID' |
| `-mic/--mut_id_column` | STR |突变 ID 列 |无 |
| `-k/--k` |整数 |侧翼序列的长度 | 30|
| `-o/--out` | STR |输出FASTA文件路径|标准输出|
| `-q/--quiet` |旗帜|抑制进步|假|

**返回：** FASTA 格式的突变序列

---

### gget gpt
使用 OpenAI 的 API 生成文本。

**设置：** 需要 `gget setup gpt` 和 OpenAI API 密钥

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `prompt` | STR |生成文本输入 |必填|
| `api_key` | STR | OpenAI API 密钥 |必填|
| `model` | STR | OpenAI 模型名称 | gpt-3.5-涡轮|
| `temperature` |浮动|采样温度（0-2）| 1.0 |
| `top_p` |浮动|细胞核取样| 1.0 |
| `max_tokens` |整数 |生成的最大代币 |无 |
| `frequency_penalty` |浮动|频率惩罚 (0-2) | 0 |
| `presence_penalty` |浮动|出场惩罚 (0-2) | 0 |

**重要提示：** 免费套餐期限为 3 个月。设置计费限额。

**返回：** 生成的文本字符串

---

### gget 设置
安装/下载模块的依赖项。

**参数：**
|参数|类型 |描述 |默认 |
|------------|------|-------------|---------|
| `module` | STR |模块名称|必填|
| `-o/--out` | STR |输出文件夹（仅限 elm）|包安装文件夹 |
| `-q/--quiet` |旗帜|抑制进步|假|

**需要设置的模块：**
- `alphafold` - 下载 ~4GB 模型参数
- `cellxgene` - 安装 cellxgene-census
- `elm` - 下载本地 ELM 数据库
- `gpt` - 配置 OpenAI 集成

**返回：** None（安装依赖项）