<!-- 此文件由机器翻译自 query_syntax.md -->

# UniProt 查询语法参考

用于构建复杂搜索的 UniProt 搜索查询语法的综合指南。

## 基本语法

### 简单查询
```
insulin
kinase
```

### 特定领域的搜索
<<<代码块_1>>>

## 布尔运算符

### AND（两个术语都必须出现）
<<<代码块_2>>>

### OR（任一术语都可以出现）
<<<代码块_3>>>

### NOT（排除术语）
<<<代码块_4>>>

### 用括号分组
<<<代码块_5>>>

## 常用搜索字段

### 身份识别
- `accession:P12345` - UniProt 登录号
- `id:INSR_HUMAN` - 条目名称
- `gene:BRCA1` - 基因名称
- `gene_exact:BRCA1` - 基因名称精确匹配

### 生物/分类学
- `organism_name:human` - 生物体名称
- `organism_name:"Homo sapiens"` - 准确的生物体名称（对多个单词使用引号）
- `organism_id:9606` - NCBI 分类 ID
- `taxonomy_id:9606` - 与organism_id相同
- `taxonomy_name:"Homo sapiens"` - 分类名称

### 蛋白质信息
- `protein_name:insulin` - 蛋白质名称
- `protein_name:"insulin receptor"` - 确切的蛋白质名称
- `reviewed:true` - 仅 Swiss-Prot（已审核）条目
- `reviewed:false` - 仅 TrEMBL（未经审核）条目

### 序列属性
- `length:[100 TO 500]` - 序列长度范围
- `mass:[50000 TO 100000]` - 道尔顿分子质量
- `sequence:MVLSPADKTNVK` - 精确序列匹配
- `fragment:false` - 排除片段序列

### 基因本体论（GO）
- `go:0005515` - GO 术语 ID（0005515 = 蛋白质结合）
- `go_f:* ` - 任何分子函数
- `go_p:*` - 任何生物过程
- `go_c:*` - 任何细胞成分

### 注释
- `annotation:(type:signal)` - 有信号肽注释
- `annotation:(type:transmem)` - 具有跨膜区域
- `cc_function:*` - 有函数注释
- `cc_interaction:*` - 有交互评论
- `ft_domain:*` - 具有域功能

### 数据库交叉引用
- `xref:pdb` - 具有 PDB 结构
- `xref:ensembl` - 有 Ensembl 参考
- `database:pdb` - 与外部参照相同
- `database:(type:pdb)` - 替代语法

### 蛋白质家族和结构域
- `family:"protein kinase"` - 蛋白质家族
- `keyword:"Protein kinase"` - 关键字注释
- `cc_similarity:*` - 有相似注释

## 范围查询

### 数字范围
<<<代码块_6>>>

### 日期范围
```
created:[2023-01-01 TO 2023-12-31]
modified:[2024-01-01 TO *]
```

## 通配符

### 单字符 (?)
```
gene:BRCA?      # Matches BRCA1, BRCA2, etc.
```

### 多个字符 (*)
```
gene:BRCA*      # Matches BRCA1, BRCA2, BRCA1P1, etc.
protein_name:kinase*
organism_name:Homo*
```

## 高级搜索

### 存在查询
```
cc_function:*              # Has any function annotation
ft_domain:*                # Has any domain feature
xref:pdb                   # Has PDB structure
```

### 组合复杂查询
```
# Human reviewed kinases with PDB structure
(protein_name:kinase OR family:kinase) AND organism_id:9606 AND reviewed:true AND xref:pdb

# Cancer-related proteins excluding mice
(disease:cancer OR keyword:cancer) NOT organism_name:mouse

# Membrane proteins with signal peptides
annotation:(type:transmem) AND annotation:(type:signal) AND reviewed:true

# Recently updated human proteins
organism_id:9606 AND modified:[2024-01-01 TO *] AND reviewed:true
```

## 特定领域的示例

### 蛋白质名称
```
protein_name:"insulin receptor"    # Exact phrase
protein_name:insulin*              # Starts with insulin
recommended_name:insulin           # Recommended name only
alternative_name:insulin           # Alternative names only
```

### 基因
```
gene:BRCA1                        # Gene symbol
gene_exact:BRCA1                  # Exact gene match
olnName:BRCA1                     # Ordered locus name
orfName:BRCA1                     # ORF name
```

### 生物体
```
organism_name:human               # Common name
organism_name:"Homo sapiens"      # Scientific name
organism_id:9606                  # Taxonomy ID
lineage:primates                  # Taxonomic lineage
```

### 特点
```
ft_signal:*                       # Signal peptide
ft_transmem:*                     # Transmembrane region
ft_domain:"Protein kinase"        # Specific domain
ft_binding:*                      # Binding site
ft_site:*                         # Any site
```

### 评论 (cc_)
```
cc_function:*                     # Function description
cc_catalytic_activity:*           # Catalytic activity
cc_pathway:*                      # Pathway involvement
cc_interaction:*                  # Protein interactions
cc_subcellular_location:*         # Subcellular location
cc_tissue_specificity:*           # Tissue specificity
cc_disease:cancer                 # Disease association
```

## 提示和最佳实践

1. **对精确短语使用引号**：`organism_name:"Homo sapiens"` 而不是 `organism_name:Homo sapiens`

2. **按审核状态过滤**：添加 `AND reviewed:true` 以获取高质量的 Swiss-Prot 条目

3. **仔细组合通配符**：`*kinase*`可能太宽泛； `kinase*` 更具体

4. **复杂逻辑使用括号**：`(A OR B) AND (C OR D)`比`A OR B AND C OR D`更清晰

5. **包含数字范围**：`length:[100 TO 500]` 包括 100 和 500

6. **字段前缀**：学习常用前缀：
   - `cc_` = 评论
   - `ft_` = 功能
   - `go_` = 基因本体
   - `xref_` = 交叉引用

7. **检查字段名称**：使用 API 的 `/configure/uniprotkb/result-fields` 端点查看所有可用字段

## 查询验证

使用以下方法测试查询：
- **网络界面**：https://www.uniprot.org/uniprotkb
- **API**：https://rest.uniprot.org/uniprotkb/search?query=YOUR_QUERY
- **API 文档**：https://www.uniprot.org/help/query-fields

## 常见模式

### 寻找特征明确的蛋白质
```
reviewed:true AND xref:pdb AND cc_function:*
```

### 寻找疾病相关蛋白质
```
cc_disease:* AND organism_id:9606 AND reviewed:true
```

### 寻找有实验证据的蛋白质
```
existence:"Evidence at protein level" AND reviewed:true
```

### 寻找分泌蛋白
```
cc_subcellular_location:secreted AND reviewed:true
```

### 寻找药物靶点
```
keyword:"Pharmaceutical" OR keyword:"Drug target"
```

## 资源

- 完整查询字段参考：https://www.uniprot.org/help/query-fields
- API查询文档：https://www.uniprot.org/help/api_queries
- 文本搜索文档：https://www.uniprot.org/help/text-search