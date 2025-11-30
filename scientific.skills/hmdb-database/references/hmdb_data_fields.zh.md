<!-- 此文件由机器翻译自 hmdb_data_fields.md -->

# HMDB 数据字段参考

本文档提供有关 HMDB 代谢物条目中可用数据字段的详细信息。

## 代谢物条目结构

每个 HMDB 代谢物条目包含 130 多个数据字段，分为几个类别：

### 化学数据字段

**鉴定：**
- `accession`：主 HMDB ID（例如 HMDB0000001）
- `secondary_accessions`：合并条目的先前 HMDB ID
- `name`：主要代谢物名称
- `synonyms`：备用名称和通用名称
- `chemical_formula`：分子式（例如，C6H12O6）
- `average_molecular_weight`：平均分子量（道尔顿）
- `monoisotopic_molecular_weight`：单同位素分子量

**结构表示：**
- `smiles`：简化的分子输入行输入系统字符串
- `inchi`：国际化学品标识符字符串
- `inchikey`：散列 InChI 用于快速查找
- `iupac_name`：IUPAC 系统名称
- `traditional_iupac`：传统的 IUPAC 名称

**化学性质：**
- `state`：物理状态（固态、液态、气态）
- `charge`：分子净电荷
- `logp`：辛醇-水分配系数（实验/预测）
- `pka_strongest_acidic`：最强酸性 pKa 值
- `pka_strongest_basic`：最强的基本 pKa 值
- `polar_surface_area`：拓扑极表面积 (TPSA)
- `refractivity`：摩尔折射率
- `polarizability`：分子极化率
- `rotatable_bond_count`：可旋转键的数量
- `acceptor_count`：氢键受体计数
- `donor_count`：氢键供体计数

**化学分类：**
- `kingdom`：化学王国（例如有机化合物）
- `super_class`：化学超类
- `class`：化学类别
- `sub_class`：化学子类
- `direct_parent`：直接化学母体
- `alternative_parents`：替代父分类
- `substituents`：存在化学取代基
- `description`：化合物的文本描述

### 生物数据字段

**代谢物来源：**
- `origin`：代谢物来源（内源性、外源性、药物代谢物、食物成分）
- `biofluid_locations`：发现的生物液体（血液、尿液、唾液、脑脊液等）
- `tissue_locations`：发现的组织（肝、肾、脑、肌肉等）
- `cellular_locations`：亚细胞位置（细胞质、线粒体、膜等）

**生物样本信息：**
- `biospecimen`：生物样本类型
- `status`：检测状态（检测到、预期、预测）
- `concentration`：带单位的浓度范围
- `concentration_references`：浓度数据的引用

**正常和异常浓度：**
对于每种生物流体（血液、尿液、唾液、脑脊液、粪便、汗液）：
- 正常浓度值和范围
- 单位（μM、mg/L 等）
- 年龄和性别考虑因素
- 浓度指标异常
- 临床意义

### 途径和酶信息

**代谢途径：**
- `pathways`：相关代谢途径列表
  - 路径名称
  - SMPDB ID（小分子通路数据库ID）
  - KEGG 通路 ID
  - 衔接类别

**酶促反应：**
- `protein_associations`：酶和转运蛋白
  - 蛋白质名称
  - 基因名称
  - Uniprot ID
  - GenBank ID
  - 蛋白质类型（酶、转运蛋白、载体等）
  - 酶反应
  - 酶动力学（Km值）

**生化背景：**
- `reactions`：涉及代谢物的生化反应
- `reaction_enzymes`：酶催化反应
- `cofactors`：必需的辅因子
- `inhibitors`：已知的酶抑制剂

### 疾病和生物标志物协会

**疾病链接：**
- `diseases`：相关疾病和状况
  - 疾病名称
  - OMIM ID（人类在线孟德尔遗传）
  - 疾病类别
  - 参考文献和证据

**生物标志物信息：**
- `biomarker_status`：化合物是否是已知的生物标志物
- `biomarker_applications`：临床应用
- `biomarker_for`：用作生物标志物的疾病或状况

### 光谱数据

**核磁共振谱：**
- `nmr_spectra`：核磁共振数据
  - 光谱类型（1D 1H、13C、2D COSY、HSQC 等）
  - 光谱仪频率（MHz）
  - 使用溶剂
  - 温度
  - 酸碱度
  - 具有化学位移和多重性的峰列表
  - FID（自由感应衰减）文件

**质谱分析：**
- `ms_spectra`：质谱数据
  - 光谱类型（MS、MS-MS、LC-MS、GC-MS）
  - 电离模式（正、负、中性）
  - 碰撞能量
  - 仪器类型
  - 峰列表（m/z、强度、注释）
  - 预测与实验标志

**色谱法：**
- `chromatography`：色谱属性
  - 保留时间
  - 柱型
  - 流动相
  - 方法详情

### 外部数据库链接

**数据库交叉引用：**
- `kegg_id`：KEGG 化合物 ID
- `pubchem_compound_id`：PubChem CID
- `pubchem_substance_id`：PubChem SID
- `chebi_id`：具有生物意义的化学实体 ID
- `chemspider_id`：ChemSpider ID
- `drugbank_id`：DrugBank 登录（如果适用）
- `foodb_id`：FooDB ID（如果是食品成分）
- `knapsack_id`：KNApSAcK ID
- `metacyc_id`：MetaCyc ID
- `bigg_id`：BiGG 型号 ID
- `wikipedia_id`：维基百科页面链接
- `metlin_id`：METLIN ID
- `vmh_id`：虚拟代谢人类 ID
- `fbonto_id`：FlyBase本体ID

**蛋白质数据库链接：**
- `uniprot_id`：相关蛋白的 UniProt 登录
- `genbank_id`：相关基因的 GenBank ID
- `pdb_id`：蛋白质结构的蛋白质数据库 ID

### 文献和证据

**参考资料：**
- `general_references`：有关代谢物的一般参考
  - 考研 ID
  - 参考文本
  - 引文
- `synthesis_reference`：合成方法和参考
- `protein_references`：蛋白质关联参考
- `pathway_references`：路径参与的参考

### 本体和分类

**本体术语：**
- `ontology_terms`：相关本体分类
  - 术语名称
  - 本体来源（ChEBI、MeSH 等）
  - 术语ID
  - 定义

### 数据质量和来源

**元数据：**
- `creation_date`：创建日期条目
- `update_date`：上次更新日期条目
- `version`：HMDB版本号
- `status`：条目状态（检测到、预期、预测）
- `evidence`：检测/存在的证据级别

## XML 结构示例

下载 XML 格式的 HMDB 数据时，结构遵循以下模式：

```xml
<metabolite>
  <accession>HMDB0000001</accession>
  <name>1-Methylhistidine</name>
  <chemical_formula>C7H11N3O2</chemical_formula>
  <average_molecular_weight>169.1811</average_molecular_weight>
  <monoisotopic_molecular_weight>169.085126436</monoisotopic_molecular_weight>
  <smiles>CN1C=NC(CC(=O)O)=C1</smiles>
  <inchi>InChI=1S/C7H11N3O2/c1-10-4-8-3-5(10)2-7(11)12/h3-4H,2H2,1H3,(H,11,12)</inchi>
  <inchikey>BRMWTNUJHUMWMS-UHFFFAOYSA-N</inchikey>

  <biospecimen_locations>
    <biospecimen>Blood</biospecimen>
    <biospecimen>Urine</biospecimen>
  </biospecimen_locations>

  <pathways>
    <pathway>
      <name>Histidine Metabolism</name>
      <smpdb_id>SMP0000044</smpdb_id>
      <kegg_map_id>map00340</kegg_map_id>
    </pathway>
  </pathways>

  <diseases>
    <disease>
      <name>Carnosinemia</name>
      <omim_id>212200</omim_id>
    </disease>
  </diseases>

  <normal_concentrations>
    <concentration>
      <biospecimen>Blood</biospecimen>
      <concentration_value>3.8</concentration_value>
      <concentration_units>uM</concentration_units>
    </concentration>
  </normal_concentrations>
</metabolite>
```

## 查询特定字段

以编程方式处理 HMDB 数据时：

**用于代谢物鉴定：**
- 按`accession`、`name`、`synonyms`、`inchi`、`smiles`查询

**对于化学相似性：**
- 使用`smiles`、`inchi`、`inchikey`、`molecular_weight`、`chemical_formula`

**对于生物标志物发现：**
- 按`diseases`、`biomarker_status`、`normal_concentrations`、`abnormal_concentrations`过滤

**对于路径分析：**
- 提取`pathways`、`protein_associations`、`reactions`

**对于光谱匹配：**
- 与 `nmr_spectra`、`ms_spectra` 峰值列表进行比较

**对于跨数据库集成：**
- 使用外部 ID 进行映射：`kegg_id`、`pubchem_compound_id`、`chebi_id` 等。

## 字段完整性

并非每种代谢物的所有字段均已填充：

- **高度完整的字段**（> 90% 的条目）：加入、名称、化学式、分子重量、微笑、英寸
- **中等完整** (50-90%)：生物样本位置、组织位置、路径
- **不同程度地完整** (10-50%)：浓度数据、疾病关联、蛋白质关联
- **稀疏完整** (<10%)：实验 NMR/MS 谱，详细的动力学数据

预测和计算数据（例如，预测的 MS 谱图、预测的浓度）补充了可用的实验数据。