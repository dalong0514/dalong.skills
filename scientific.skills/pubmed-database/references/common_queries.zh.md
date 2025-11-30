<!-- 此文件由机器翻译自 common_queries.md -->

# 常见的 PubMed 查询模式

本参考文献提供了针对各种研究场景的常见 PubMed 搜索模式的实际示例。

## 一般研究查询

### 查找某个主题的最新研究
```
breast cancer[tiab] AND 2023:2024[dp]
```

### 对某个主题的系统评论
<<<代码块_1>>>

### 荟萃分析
<<<代码块_2>>>

### 临床试验
<<<代码块_3>>>

### 寻找指南
<<<代码块_4>>>

## 疾病特定查询

### 癌症研究
<<<代码块_5>>>

### 心血管疾病
<<<代码块_6>>>

### 传染病
```
# COVID-19 research
COVID-19[tiab] AND (vaccine[tiab] OR vaccination[tiab]) AND 2023:2024[dp]

# Antibiotic resistance
(antibiotic resistance[tiab] OR drug resistance, bacterial[mh]) AND systematic review[pt]

# Tuberculosis treatment
tuberculosis[mh]/drug therapy AND (multidrug-resistant[tiab] OR MDR-TB[tiab])
```

### 神经系统疾病
```
# Alzheimer's disease
alzheimer disease[mh] AND (diagnosis[sh] OR biomarkers[tiab]) AND 2020:2024[dp]

# Parkinson's disease treatment
parkinson disease[mh] AND treatment[tiab] AND clinical trial[pt]

# Multiple sclerosis
multiple sclerosis[mh] AND disease modifying[tiab] AND review[pt]
```

### 糖尿病
```
# Type 2 diabetes management
diabetes mellitus, type 2[mh] AND (lifestyle[tiab] OR diet[tiab]) AND randomized controlled trial[pt]

# Diabetes complications
diabetes mellitus[mh] AND (complications[sh] OR diabetic neuropathy[mh])

# New diabetes drugs
diabetes mellitus, type 2[mh] AND (GLP-1[tiab] OR SGLT2[tiab]) AND 2022:2024[dp]
```

## 药物和治疗研究

### 药物功效研究
```
# Compare two drugs
(drug A[nm] OR drug B[nm]) AND condition[mh] AND comparative effectiveness[tiab]

# Drug side effects
medication name[nm] AND (adverse effects[sh] OR side effects[tiab])

# Drug combination therapy
(aspirin[nm] AND clopidogrel[nm]) AND acute coronary syndrome[mh]
```

### 治疗比较
```
# Surgery vs medication
condition[mh] AND (surgery[tiab] OR surgical[tiab]) AND (medication[tiab] OR drug therapy[sh]) AND comparative study[pt]

# Different surgical approaches
procedure[tiab] AND (laparoscopic[tiab] OR open surgery[tiab]) AND outcomes[tiab]
```

### 替代医学
```
# Herbal supplements
(herbal medicine[mh] OR phytotherapy[mh]) AND condition[tiab] AND clinical trial[pt]

# Acupuncture
acupuncture[mh] AND pain[tiab] AND randomized controlled trial[pt]
```

## 诊断研究

### 诊断测试
```
# Sensitivity and specificity
test name[tiab] AND condition[tiab] AND (sensitivity[tiab] AND specificity[tiab])

# Diagnostic imaging
(MRI[tiab] OR magnetic resonance imaging[tiab]) AND brain tumor[tiab] AND diagnosis[sh]

# Lab test evaluation
biomarker name[tiab] AND disease[tiab] AND (diagnostic[tiab] OR screening[tiab])
```

### 筛选计划
```
# Cancer screening
cancer type[tiab] AND screening[tiab] AND (cost effectiveness[tiab] OR benefit[tiab])

# Population screening
condition[tiab] AND mass screening[mh] AND public health[tiab]
```

## 针对特定人群的查询

### 儿科研究
```
# Children with specific condition
condition[tiab] AND (child[mh] OR pediatric[tiab]) AND treatment[tiab]

# Age-specific
disease[tiab] AND (infant[mh] OR child, preschool[mh])

# Pediatric dosing
drug name[nm] AND pediatric[tiab] AND (dosing[tiab] OR dose[tiab])
```

### 老年医学研究
```
# Elderly population
condition[tiab] AND (aged[mh] OR elderly[tiab] OR geriatric[tiab])

# Aging and disease
aging[mh] AND disease[tiab] AND mechanism[tiab]

# Polypharmacy
polypharmacy[tiab] AND elderly[tiab] AND adverse effects[tiab]
```

### 孕妇
```
# Pregnancy and medications
drug name[nm] AND (pregnancy[mh] OR pregnant women[tiab]) AND safety[tiab]

# Pregnancy complications
pregnancy complication[tiab] AND management[tiab]
```

### 针对性别的研究
```
# Female-specific
condition[tiab] AND female[mh] AND hormones[tiab]

# Male-specific
disease[tiab] AND male[mh] AND risk factors[tiab]

# Sex differences
condition[tiab] AND (sex factors[mh] OR gender differences[tiab])
```

## 流行病学和公共卫生

### 患病率研究
```
disease[tiab] AND (prevalence[tiab] OR epidemiology[sh]) AND country/region[tiab]
```

### 发生率研究
```
condition[tiab] AND incidence[tiab] AND population[tiab] AND 2020:2024[dp]
```

### 风险因素
```
disease[mh] AND (risk factors[mh] OR etiology[sh]) AND cohort study[tiab]
```

### 全球健康
```
disease[tiab] AND (developing countries[mh] OR low income[tiab]) AND burden[tiab]
```

### 健康差异
```
condition[tiab] AND (health disparities[tiab] OR health equity[tiab]) AND minority groups[tiab]
```

## 方法特定的查询

### 研究方法

#### 队列研究
```
condition[tiab] AND cohort study[tiab] AND prospective[tiab]
```

#### 病例对照研究
```
disease[tiab] AND case-control studies[mh] AND risk factors[tiab]
```

#### 横断面研究
```
condition[tiab] AND cross-sectional studies[mh] AND prevalence[tiab]
```

### 统计方法
```
# Machine learning in medicine
(machine learning[tiab] OR artificial intelligence[tiab]) AND diagnosis[tiab] AND validation[tiab]

# Bayesian analysis
condition[tiab] AND bayes theorem[mh] AND clinical decision[tiab]
```

### 遗传和分子研究
```
# GWAS studies
disease[tiab] AND (genome-wide association study[tiab] OR GWAS[tiab])

# Gene expression
gene name[tiab] AND (gene expression[mh] OR mRNA[tiab]) AND disease[tiab]

# Proteomics
condition[tiab] AND proteomics[mh] AND biomarkers[tiab]

# CRISPR research
CRISPR[tiab] AND (gene editing[tiab] OR genome editing[tiab]) AND 2020:2024[dp]
```

## 作者和机构查询

### 查找特定作者的作品
```
# Single author
smith ja[au] AND cancer[tiab] AND 2023:2024[dp]

# First author only
jones m[1au] AND cardiology[tiab]

# Multiple authors from same group
(smith ja[au] OR jones m[au] OR wilson k[au]) AND research topic[tiab]
```

### 针对机构的研究
```
# University affiliation
harvard[affil] AND cancer research[tiab] AND 2023:2024[dp]

# Hospital research
"mayo clinic"[affil] AND clinical trial[pt]

# Country-specific
japan[affil] AND robotics[tiab] AND surgery[tiab]
```

## 期刊特定查询

### 高影响力期刊
```
# Specific journal
nature[ta] AND genetics[tiab] AND 2024[dp]

# Multiple journals
(nature[ta] OR science[ta] OR cell[ta]) AND immunology[tiab]

# Journal with ISSN
0028-4793[issn] AND clinical trial[pt]
```

## 引文和参考文献查询

### 查找特定文章
```
# By PMID
12345678[pmid]

# By DOI
10.1056/NEJMoa123456[doi]

# By first author and year
smith ja[1au] AND 2023[dp] AND cancer[tiab]
```

### 查找被引著作
```
# Related articles
Similar Articles feature from any PubMed result

# By keyword in references
Use "Cited by" links when available
```

## 高级组合查询

### 综合文献综述
```
(disease name[tiab] OR disease name[mh]) AND
((treatment[tiab] OR therapy[tiab] OR management[tiab]) OR
(diagnosis[tiab] OR screening[tiab]) OR
(epidemiology[tiab] OR prevalence[tiab])) AND
(systematic review[pt] OR meta-analysis[pt] OR review[pt]) AND
2019:2024[dp] AND english[la]
```

### 精准医学查询
```
(precision medicine[tiab] OR personalized medicine[tiab] OR pharmacogenomics[mh]) AND
cancer[tiab] AND
(biomarkers[tiab] OR genetic testing[tiab]) AND
clinical application[tiab] AND
2020:2024[dp]
```

### 转化研究
```
(basic science[tiab] OR bench to bedside[tiab] OR translational medical research[mh]) AND
disease[tiab] AND
(clinical trial[pt] OR clinical application[tiab]) AND
2020:2024[dp]
```

## 质量过滤器

### 高质量证据
```
condition[tiab] AND
(randomized controlled trial[pt] OR systematic review[pt] OR meta-analysis[pt]) AND
humans[mh] AND
english[la] AND
2020:2024[dp]
```

### 免费全文文章
```
topic[tiab] AND free full text[sb] AND 2023:2024[dp]
```

### 带摘要的文章
```
condition[tiab] AND hasabstract[text] AND review[pt]
```

## 保持最新状态

### 最新出版物
```
topic[tiab] AND 2024[dp] AND english[la]
```

### 预印本和抢先体验
```
topic[tiab] AND (epub ahead of print[tiab] OR publisher[sb])
```

### 设置警报
```
# Create search and save to My NCBI
# Enable email alerts for new matching articles
topic[tiab] AND (randomized controlled trial[pt] OR systematic review[pt])
```

## COVID-19 特定查询

### 疫苗研究
```
(COVID-19[tiab] OR SARS-CoV-2[tiab]) AND
(vaccine[tiab] OR vaccination[tiab]) AND
(efficacy[tiab] OR effectiveness[tiab]) AND
2023:2024[dp]
```

### 长新冠病毒
```
(long covid[tiab] OR post-acute covid[tiab] OR PASC[tiab]) AND
(symptoms[tiab] OR treatment[tiab])
```

### 新冠肺炎治疗
```
COVID-19[tiab] AND
(antiviral[tiab] OR monoclonal antibody[tiab] OR treatment[tiab]) AND
randomized controlled trial[pt]
```

## 构建查询的技巧

### 1.PICO框架
使用 PICO（群体、干预、比较、结果）构建临床查询：

```
P: diabetes mellitus, type 2[mh]
I: metformin[nm]
C: lifestyle modification[tiab]
O: glycemic control[tiab]

Query: diabetes mellitus, type 2[mh] AND (metformin[nm] OR lifestyle modification[tiab]) AND glycemic control[tiab]
```

### 2.迭代细化
从广泛开始，审查结果，完善：
```
1. diabetes → too broad
2. diabetes mellitus type 2 → better
3. diabetes mellitus, type 2[mh] AND metformin[nm] → more specific
4. diabetes mellitus, type 2[mh] AND metformin[nm] AND randomized controlled trial[pt] → focused
```

### 3. 使用搜索历史记录
在高级搜索中合并以前的搜索：
```
#1: diabetes mellitus, type 2[mh]
#2: cardiovascular disease[mh]
#3: #1 AND #2 AND risk factors[tiab]
```

### 4. 保存有效搜索
创建我的 NCBI 帐户以保存成功的查询以供将来使用并设置自动警报。