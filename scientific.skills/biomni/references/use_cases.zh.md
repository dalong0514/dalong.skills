<!-- 此文件由机器翻译自 use_cases.md -->

# Biomni 用例和示例

展示跨生物医学研究领域的 biomni 的综合示例。

## 目录

1. [CRISPR筛选和基因编辑](#crispr-screening-and-gene-editing)
2. [单细胞RNA-seq分析](#single-cell-rna-seq-analysis)
3. [药物发现和 ADMET](#drug-discovery-and-admet)
4. [GWAS 和遗传分析](#gwas-and-genic-analysis)
5. [临床基因组学和诊断](#clinical-genomics-and-diagnostics)
6. [蛋白质结构与功能](#蛋白质结构与功能)
7.【文献与知识综合】(#文献与知识综合)
8. [多组学整合](#multi-omics-integration)

---

## CRISPR 筛选和基因编辑

### 示例 1：全基因组 CRISPR 筛选设计

**任务：** 设计 CRISPR 敲除筛选来识别调节自噬的基因。

```python
from biomni.agent import A1

agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Design a genome-wide CRISPR knockout screen to identify genes regulating
autophagy in HEK293 cells.

Requirements:
1. Generate comprehensive sgRNA library targeting all protein-coding genes
2. Design 4 sgRNAs per gene with optimal on-target and minimal off-target scores
3. Include positive controls (known autophagy regulators: ATG5, BECN1, ULK1)
4. Include negative controls (non-targeting sgRNAs)
5. Prioritize genes based on:
   - Existing autophagy pathway annotations
   - Protein-protein interactions with known autophagy factors
   - Expression levels in HEK293 cells
6. Output sgRNA sequences, scores, and gene prioritization rankings

Provide analysis as Python code and interpret results.
""")

agent.save_conversation_history("autophagy_screen_design.pdf")
```

**预期输出：**
- sgRNA 文库，包含约 80,000 个向导（每个基因 4 个 × 约 20,000 个基因）
- 每个 sgRNA 的中靶和脱靶分数
- 基于通路富集的优先基因列表
- 文库设计的质量控制指标

### 示例 2：CRISPR 脱靶预测

<<<代码块_1>>>

### 示例 3：屏幕点击分析

<<<代码块_2>>>

---

## 单细胞 RNA-seq 分析

### 示例 1：单元格类型注释

**任务：** 分析单细胞 RNA-seq 数据并注释细胞群。

<<<代码块_3>>>

### 示例 2：差异表达分析

<<<代码块_4>>>

### 示例 3：轨迹分析

<<<代码块_5>>>

---

## 药物发现和 ADMET

### 示例 1：ADMET 属性预测

**任务：** 预测候选药物的 ADMET 特性。

<<<代码块_6>>>

### 示例 2：目标识别

```python
result = agent.go("""
Identify potential protein targets for Alzheimer's disease drug development.

Tasks:
1. Query GWAS data for Alzheimer's-associated genes
2. Identify genes with druggable domains (kinases, GPCRs, ion channels, etc.)
3. Check for brain expression patterns
4. Assess disease relevance via literature mining
5. Evaluate existing chemical probe availability
6. Rank targets by:
   - Genetic evidence strength
   - Druggability
   - Lack of existing therapies
7. Suggest target validation experiments
""")
```

### 示例 3：虚拟筛选

```python
result = agent.go("""
Perform virtual screening for EGFR kinase inhibitors.

Database: ZINC15 lead-like subset (~6M compounds)
Target: EGFR kinase domain (PDB: 1M17)

Workflow:
1. Prepare protein structure (remove waters, add hydrogens)
2. Define binding pocket (based on erlotinib binding site)
3. Generate pharmacophore model from known EGFR inhibitors
4. Filter ZINC database by:
   - Molecular weight: 200-500 Da
   - LogP: 0-5
   - Lipinski's rule of five
   - Pharmacophore match
5. Dock top 10,000 compounds
6. Score by docking energy and predicted binding affinity
7. Select top 100 for further analysis
8. Predict ADMET properties for top hits
9. Recommend top 10 compounds for experimental validation
""")
```

---

## GWAS 和遗传分析

### 示例1：GWAS汇总统计分析

**任务：** 解释 GWAS 结果并识别因果基因。

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Analyze GWAS summary statistics for Type 2 Diabetes.

Input file: t2d_gwas_summary.txt
Columns: CHR, BP, SNP, P, OR, BETA, SE, A1, A2

Analysis steps:
1. Identify genome-wide significant variants (P < 5e-8)
2. Perform LD clumping to identify independent signals
3. Map variants to genes using:
   - Nearest gene
   - eQTL databases (GTEx)
   - Hi-C chromatin interactions
4. Prioritize causal genes using multiple evidence:
   - Fine-mapping scores
   - Coding variant consequences
   - Gene expression in relevant tissues (pancreas, liver, adipose)
   - Pathway enrichment
5. Identify druggable targets among causal genes
6. Compare with known T2D genes and highlight novel associations
7. Generate Manhattan plot, QQ plot, and gene prioritization table
""")

agent.save_conversation_history("t2d_gwas_analysis.pdf")
```

### 示例 2：多基因风险评分

```python
result = agent.go("""
Develop and validate polygenic risk score (PRS) for coronary artery disease (CAD).

Training GWAS: CAD_discovery_summary_stats.txt (N=180,000)
Validation cohort: CAD_validation_genotypes.vcf (N=50,000)

Tasks:
1. Select variants for PRS using p-value thresholding (P < 1e-5)
2. Perform LD clumping (r² < 0.1, 500kb window)
3. Calculate PRS weights from GWAS betas
4. Compute PRS for validation cohort individuals
5. Evaluate PRS performance:
   - AUC for CAD case/control discrimination
   - Odds ratios across PRS deciles
   - Compare to traditional risk factors (age, sex, BMI, smoking)
6. Assess PRS calibration and create risk stratification plot
7. Identify high-risk individuals (top 5% PRS)
""")
```

### 示例 3：变异致病性预测

```python
result = agent.go("""
Predict pathogenicity of rare coding variants in candidate disease genes.

Variants (VCF format):
- chr17:41234451:A>G (BRCA1 p.Arg1347Gly)
- chr2:179428448:C>T (TTN p.Trp13579*)
- chr7:117188679:G>A (CFTR p.Gly542Ser)

For each variant, assess:
1. In silico predictions (SIFT, PolyPhen2, CADD, REVEL)
2. Population frequency (gnomAD)
3. Evolutionary conservation (PhyloP, PhastCons)
4. Protein structure impact (using AlphaFold structures)
5. Functional domain location
6. ClinVar annotations (if available)
7. Literature evidence
8. ACMG/AMP classification criteria

Provide pathogenicity classification (benign, likely benign, VUS, likely pathogenic, pathogenic) with supporting evidence.
""")
```

---

## 临床基因组学和诊断

### 示例1：罕见病诊断

**任务：** 通过全外显子组测序诊断罕见遗传病。

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Analyze whole exome sequencing (WES) data for rare disease diagnosis.

Patient phenotypes (HPO terms):
- HP:0001250 (Seizures)
- HP:0001249 (Intellectual disability)
- HP:0001263 (Global developmental delay)
- HP:0001252 (Hypotonia)

VCF file: patient_trio.vcf (proband + parents)

Analysis workflow:
1. Variant filtering:
   - Quality filters (QUAL > 30, DP > 10, GQ > 20)
   - Frequency filters (gnomAD AF < 0.01)
   - Functional impact (missense, nonsense, frameshift, splice site)

2. Inheritance pattern analysis:
   - De novo variants
   - Autosomal recessive (compound het, homozygous)
   - X-linked

3. Phenotype-driven prioritization:
   - Match candidate genes to HPO terms
   - Use HPO-gene associations
   - Check gene expression in relevant tissues (brain)

4. Variant pathogenicity assessment:
   - In silico predictions
   - ACMG classification
   - Literature evidence

5. Generate diagnostic report with:
   - Top candidate variants
   - Supporting evidence
   - Functional validation suggestions
   - Genetic counseling recommendations
""")

agent.save_conversation_history("rare_disease_diagnosis.pdf")
```

### 示例 2：癌症基因组分析

```python
result = agent.go("""
Analyze tumor-normal paired sequencing for cancer genomics.

Files:
- tumor_sample.vcf (somatic variants)
- tumor_rnaseq.bam (gene expression)
- tumor_cnv.seg (copy number variants)

Analysis:
1. Identify driver mutations:
   - Known cancer genes (COSMIC, OncoKB)
   - Recurrent hotspot mutations
   - Truncating mutations in tumor suppressors

2. Analyze mutational signatures:
   - Decompose signatures (COSMIC signatures)
   - Identify mutagenic processes

3. Copy number analysis:
   - Identify amplifications and deletions
   - Focal vs. arm-level events
   - Assess oncogene amplifications and TSG deletions

4. Gene expression analysis:
   - Identify outlier gene expression
   - Fusion transcript detection
   - Pathway dysregulation

5. Therapeutic implications:
   - Match alterations to FDA-approved therapies
   - Identify clinical trial opportunities
   - Predict response to targeted therapies

6. Generate precision oncology report
""")
```

### 示例 3：药物基因组学

```python
result = agent.go("""
Generate pharmacogenomics report for patient genotype data.

VCF file: patient_pgx.vcf

Analyze variants affecting drug metabolism:

**CYP450 genes:**
- CYP2D6 (affects ~25% of drugs)
- CYP2C19 (clopidogrel, PPIs, antidepressants)
- CYP2C9 (warfarin, NSAIDs)
- CYP3A5 (tacrolimus, immunosuppressants)

**Drug transporter genes:**
- SLCO1B1 (statin myopathy risk)
- ABCB1 (P-glycoprotein)

**Drug targets:**
- VKORC1 (warfarin dosing)
- DPYD (fluoropyrimidine toxicity)
- TPMT (thiopurine toxicity)

For each gene:
1. Determine diplotype (*1/*1, *1/*2, etc.)
2. Assign metabolizer phenotype (PM, IM, NM, RM, UM)
3. Provide dosing recommendations using CPIC/PharmGKB guidelines
4. Flag high-risk drug-gene interactions
5. Suggest alternative medications if needed

Generate patient-friendly report with actionable recommendations.
""")
```

---

## 蛋白质结构和功能

### 示例 1：AlphaFold 结构分析

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Analyze AlphaFold structure prediction for novel protein.

Protein: Hypothetical protein ABC123 (UniProt: Q9XYZ1)

Tasks:
1. Retrieve AlphaFold structure from database
2. Assess prediction quality:
   - pLDDT scores per residue
   - Identify high-confidence regions (pLDDT > 90)
   - Flag low-confidence regions (pLDDT < 50)

3. Structural analysis:
   - Identify domains using structural alignment
   - Predict fold family
   - Identify secondary structure elements

4. Functional prediction:
   - Search for structural homologs in PDB
   - Identify conserved functional sites
   - Predict binding pockets
   - Suggest possible ligands/substrates

5. Variant impact analysis:
   - Map disease-associated variants to structure
   - Predict structural consequences
   - Identify variants affecting binding sites

6. Generate PyMOL visualization scripts highlighting key features
""")

agent.save_conversation_history("alphafold_analysis.pdf")
```

### 示例 2：蛋白质-蛋白质相互作用预测

```python
result = agent.go("""
Predict and analyze protein-protein interactions for autophagy pathway.

Query proteins: ATG5, ATG12, ATG16L1

Analysis:
1. Retrieve known interactions from:
   - STRING database
   - BioGRID
   - IntAct
   - Literature mining

2. Predict novel interactions using:
   - Structural modeling (AlphaFold-Multimer)
   - Coexpression analysis
   - Phylogenetic profiling

3. Analyze interaction interfaces:
   - Identify binding residues
   - Assess interface properties (area, hydrophobicity)
   - Predict binding affinity

4. Functional analysis:
   - Map interactions to autophagy pathway steps
   - Identify regulatory interactions
   - Predict complex stoichiometry

5. Therapeutic implications:
   - Identify druggable interfaces
   - Suggest peptide inhibitors
   - Design disruption strategies

Generate network visualization and interaction details.
""")
```

---

## 文献与知识综合

### 示例1：系统文献综述

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Perform systematic literature review on CRISPR base editing applications.

Search query: "CRISPR base editing" OR "base editor" OR "CBE" OR "ABE"
Date range: 2016-2025

Tasks:
1. Search PubMed and retrieve relevant abstracts
2. Filter for original research articles
3. Extract key information:
   - Base editor type (CBE, ABE, dual)
   - Target organism/cell type
   - Application (disease model, therapy, crop improvement)
   - Editing efficiency
   - Off-target assessment

4. Categorize applications:
   - Therapeutic applications (by disease)
   - Agricultural applications
   - Basic research

5. Analyze trends:
   - Publications over time
   - Most studied diseases
   - Evolution of base editor technology

6. Synthesize findings:
   - Clinical trial status
   - Remaining challenges
   - Future directions

Generate comprehensive review document with citation statistics.
""")

agent.save_conversation_history("crispr_base_editing_review.pdf")
```

### 示例2：基因功能合成

```python
result = agent.go("""
Synthesize knowledge about gene function from multiple sources.

Target gene: PARK7 (DJ-1)

Integrate information from:
1. **Genetic databases:**
   - NCBI Gene
   - UniProt
   - OMIM

2. **Expression data:**
   - GTEx tissue expression
   - Human Protein Atlas
   - Single-cell expression atlases

3. **Functional data:**
   - GO annotations
   - KEGG pathways
   - Reactome
   - Protein interactions (STRING)

4. **Disease associations:**
   - ClinVar variants
   - GWAS catalog
   - Disease databases (DisGeNET)

5. **Literature:**
   - PubMed abstracts
   - Key mechanistic studies
   - Review articles

Synthesize into comprehensive gene report:
- Molecular function
- Biological processes
- Cellular localization
- Tissue distribution
- Disease associations
- Known drug targets/inhibitors
- Unresolved questions

Generate structured summary suitable for research planning.
""")
```

---

## 多组学整合

### 示例 1：多组学疾病分析

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Integrate multi-omics data to understand disease mechanism.

Disease: Alzheimer's disease
Data types:
- Genomics: GWAS summary statistics (gwas_ad.txt)
- Transcriptomics: Brain RNA-seq (controls vs AD, rnaseq_data.csv)
- Proteomics: CSF proteomics (proteomics_csf.csv)
- Metabolomics: Plasma metabolomics (metabolomics_plasma.csv)
- Epigenomics: Brain methylation array (methylation_data.csv)

Integration workflow:
1. Analyze each omics layer independently:
   - Identify significantly altered features
   - Perform pathway enrichment

2. Cross-omics correlation:
   - Correlate gene expression with protein levels
   - Link genetic variants to expression (eQTL)
   - Associate methylation with gene expression
   - Connect proteins to metabolites

3. Network analysis:
   - Build multi-omics network
   - Identify key hub genes/proteins
   - Detect disease modules

4. Causal inference:
   - Prioritize drivers vs. consequences
   - Identify therapeutic targets
   - Predict drug mechanisms

5. Generate integrative model of AD pathogenesis

Provide visualization and therapeutic target recommendations.
""")

agent.save_conversation_history("ad_multiomics_analysis.pdf")
```

### 示例 2：系统生物学建模

```python
result = agent.go("""
Build systems biology model of metabolic pathway.

Pathway: Glycolysis
Data sources:
- Enzyme kinetics (BRENDA database)
- Metabolite concentrations (literature)
- Gene expression (tissue-specific, GTEx)
- Flux measurements (C13 labeling studies)

Modeling tasks:
1. Construct pathway model:
   - Define reactions and stoichiometry
   - Parameterize enzyme kinetics (Km, Vmax, Ki)
   - Set initial metabolite concentrations

2. Simulate pathway dynamics:
   - Steady-state analysis
   - Time-course simulations
   - Sensitivity analysis

3. Constraint-based modeling:
   - Flux balance analysis (FBA)
   - Identify bottleneck reactions
   - Predict metabolic engineering strategies

4. Integrate with gene expression:
   - Tissue-specific model predictions
   - Disease vs. normal comparisons

5. Therapeutic predictions:
   - Enzyme inhibition effects
   - Metabolic rescue strategies
   - Drug target identification

Generate model in SBML format and simulation results.
""")
```

---

## 任务制定的最佳实践

### 1.具体而详细

**差：**
```python
agent.go("Analyze this RNA-seq data")
```

**好：**
```python
agent.go("""
Analyze bulk RNA-seq data from cancer vs. normal samples.

Files: cancer_rnaseq.csv (TPM values, 50 cancer, 50 normal)

Tasks:
1. Differential expression (DESeq2, padj < 0.05, |log2FC| > 1)
2. Pathway enrichment (KEGG, Reactome)
3. Generate volcano plot and top DE gene heatmap
""")
```

### 2. 包含文件路径和格式

始终指定：
- 确切的文件路径
- 文件格式（VCF、BAM、CSV、H5AD 等）
- 数据结构（列、样本 ID）

### 3. 设定明确的成功标准

定义阈值和截止值：
- 统计显着性（P < 0.05，FDR < 0.1）
- 倍数变化阈值
- 优质过滤器
- 预期产出

### 4. 请求可视化

明确要求绘图：
- 火山图、MA 图
- 热图、PCA 图
- 网络图
- 曼哈顿地块

### 5. 指定生物学背景

包括：
- 生物体（人类、小鼠等）
- 组织/细胞类型
- 疾病/状况
- 治疗详情

### 6. 请求解释

要求代理：
- 解释生物学意义
- 建议后续实验
- 确定限制
- 提供文献背景

---

## 常见模式

### 数据质量控制

```python
"""
Before analysis, perform quality control:
1. Check for missing values
2. Assess data distributions
3. Identify outliers
4. Generate QC report
Only proceed with analysis if data passes QC.
"""
```

### 迭代细化

```python
"""
Perform analysis in stages:
1. Initial exploratory analysis
2. Based on results, refine parameters
3. Focus on interesting findings
4. Generate final report

Show intermediate results for each stage.
"""
```

### 再现性

```python
"""
Ensure reproducibility:
1. Set random seeds where applicable
2. Log all parameters used
3. Save intermediate files
4. Export environment info (package versions)
5. Generate methods section for paper
"""
```

这些例子展示了 biomni 可以处理的生物医学任务的广度。根据您的具体研究问题调整模式，并始终包含足够的详细信息，以便代理自主执行。