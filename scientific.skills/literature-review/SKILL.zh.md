<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：文献综述
描述：使用多个学术数据库（PubMed、arXiv、bioRxiv、Semantic Scholar 等）进行全面、系统的文献综述。在跨生物医学、科学和技术领域进行系统文献综述、荟萃分析、研究综合或综合文献检索时，应使用此技能。创建专业格式的 Markdown 文档和 PDF，并以多种引文样式（APA、Nature、Vancouver 等）验证引文。
---

#文献综述

## 概述

按照严格的学术方法进行系统、全面的文献综述。搜索多个文献数据库，按主题综合研究结果，验证所有引用的准确性，并生成 Markdown 和 PDF 格式的专业输出文档。

该技能与数据库访问的多种科学技能（gget、生物服务、datacommons-client）集成，并提供用于引文验证、结果聚合和文档生成的专用工具。

## 何时使用此技能

在以下情况下使用此技能：
- 为研究或出版进行系统的文献综述
- 综合多个来源的特定主题的当前知识
- 进行荟萃分析或范围审查
- 撰写研究论文或论文的文献综述部分
- 调查研究领域的最新技术水平
- 确定研究差距和未来方向
- 需要经过验证的引文和专业格式

## 核心工作流程

文献综述遵循结构化、多阶段的工作流程：

### 第一阶段：规划和范围界定

1. **定义研究问题**：使用 PICO 框架（群体、干预、比较、结果）进行临床/生物医学审查
   - 示例：“与标准治疗 (C) 相比，CRISPR-Cas9 (I) 治疗镰状细胞病 (P) 的功效如何？”

2. **确定范围和目标**：
   - 定义清晰、具体的研究问题
   - 确定评论类型（叙述性、系统性、范围界定、荟萃分析）
   - 设定界限（时间段、地理范围、研究类型）

3. **制定搜索策略**：
   - 从研究问题中识别 2-4 个主要概念
   - 列出每个概念的同义词、缩写和相关术语
   - 计划布尔运算符（AND、OR、NOT）来组合术语
   - 选择至少 3 个互补数据库

4. **设置包含/排除标准**：
   - 日期范围（例如，过去 10 年：2015-2024）
   - 语言（通常为英语，或指定多语言）
   - 出版物类型（同行评审、预印本、评论）
   - 研究设计（随机对照试验、观察性试验、体外试验等）
   - 清楚地记录所有标准

### 第二阶段：系统文献检索

1. **多数据库搜索**：

   选择适合该域的数据库：

   **生物医学与生命科学：**
   - 使用 `gget` 技能：`gget search pubmed "search terms"` 用于 PubMed/PMC
   - 使用`gget`技能：`gget search biorxiv "search terms"`进行预印本
   - 对 ChEMBL、KEGG、UniProt 等使用 `bioservices` 技能。

   **一般科学文献：**
   - 通过直接 API 搜索 arXiv（物理、数学、CS、q-bio 预印本）
   - 通过 API 搜索 Semantic Scholar（2 亿多篇论文，跨学科）
   - 使用 Google Scholar 进行全面覆盖（手动或仔细抓取）

   **专业数据库：**
   - 对蛋白质结构使用 `gget alphafold`
   - 使用 `gget cosmic` 进行癌症基因组学
   - 使用 `datacommons-client` 获取人口统计/统计数据
   - 使用适合该领域的专用数据库

2. **文档搜索参数**：
   ```markdown
   ## Search Strategy

   ### Database: PubMed
   - **Date searched**: 2024-10-25
   - **Date range**: 2015-01-01 to 2024-10-25
   - **Search string**:
     ```
     （“CRISPR”[标题] 或“Cas9”[标题]）
     AND（“镰状细胞”[MeSH] 或“SCD”[标题/摘要]）
     和 2015:2024[发布日期]
     <<<代码块_1>>>

   对搜索的每个数据库重复此操作。

3. **导出并汇总结果**：
   - 从每个数据库导出 JSON 格式的结果
   - 将所有结果合并到一个文件中
   - 使用 `scripts/search_databases.py` 进行后处理：
     <<<代码块_2>>>

### 第三阶段：筛选和选择

1. **重复数据删除**：
   <<<代码块_3>>>
   - 按 DOI（主要）或标题（后备）删除重复项
   - 删除重复文档的数量

2. **标题筛选**：
   - 根据纳入/排除标准审查所有标题
   - 排除明显不相关的研究
   - 现阶段排除的文件编号

3. **摘要筛选**：
   - 阅读剩余研究的摘要
- 严格应用纳入/排除标准
   - 记录排除的原因

4. **全文筛选**：
   - 获取剩余研究的全文
   - 根据所有标准进行详细审查
   - 记录排除的具体原因
   - 记录纳入研究的最终数量

5. **创建 PRISMA 流程图**：
   <<<代码块_4>>>

### 第 4 阶段：数据提取和质量评估

1. **从每项纳入的研究中提取关键数据**：
   - 研究元数据（作者、年份、期刊、DOI）
   - 研究设计和方法
   - 样本量和总体特征
   - 主要发现和结果
   - 作者指出的局限性
   - 资金来源和利益冲突

2. **评估学习质量**：
   - **对于 RCT**：使用 Cochrane 偏差风险工具
   - **对于观察性研究**：使用纽卡斯尔-渥太华量表
   - **对于系统评价**：使用 AMSTAR 2
   - 对每项研究进行评分：高、中、低或极低质量
   - 考虑排除质量非常低的研究

3. **按主题组织**：
   - 确定研究中的 3-5 个主要主题
   - 按主题进行分组研究（研究可能出现在多个主题中）
   - 注意模式、共识和争议

### 第五阶段：综合与分析

1. **从模板创建审核文档**：
   <<<代码块_5>>>

2. **撰写主题综合**（不是逐个研究的总结）：
   - 按主题或研究问题组织结果部分
   - 综合每个主题内多项研究的结果
   - 比较和对比不同的方法和结果
   - 确定共识领域和争议点
   - 突出最有力的证据

   结构示例：
   <<<代码块_6>>>

3. **批判性分析**：
   - 评估各项研究的方法学优势和局限性
   - 评估证据的质量和一致性
   - 确定知识差距和方法差距
   - 注意未来需要研究的领域

4. **撰写讨论**：
   - 在更广泛的背景下解释研究结果
   - 讨论临床、实践或研究意义
   - 承认审查本身的局限性
   - 与之前的评论进行比较（如果适用）
   - 提出未来具体的研究方向

### 第 6 阶段：引文验证

**重要**：在最终提交之前，所有引文必须经过准确性验证。

1. **验证所有 DOI**：
   ```bash
   python scripts/verify_citations.py my_literature_review.md
   ```

   这个脚本：
   - 从文档中提取所有 DOI
   - 验证每个 DOI 正确解析
   - 从 CrossRef 检索元数据
   - 生成验证报告
   - 输出格式正确的引文

2. **审查验证报告**：
   - 检查是否有任何失败的 DOI
   - 验证作者姓名、标题和出版物详细信息是否匹配
   - 更正原始文档中的任何错误
   - 重新运行验证，直到所有引用通过

3. **引文格式保持一致**：
   - 选择一种引用样式并在全文中使用（参见`references/citation_styles.md`）
   - 常见样式：APA、Nature、Vancouver、Chicago、IEEE
   - 使用验证脚本输出正确格式化引文
   - 确保文本引用符合参考列表格式

### 第 7 阶段：文档生成

1. **生成PDF**：
   ```bash
   python scripts/generate_pdf.py my_literature_review.md \
     --citation-style apa \
     --output my_review.pdf
   ```

   选项：
   - `--citation-style`：apa、自然、芝加哥、温哥华、ieee
   - `--no-toc`：禁用目录
   - `--no-numbers`：禁用节编号
   - `--check-deps`：检查是否安装了 pandoc/xelatex

2. **查看最终输出**：
   - 检查 PDF 格式和布局
   - 验证所有部分都存在
   - 确保引文正确呈现
   - 检查数字/表格是否正确显示
   - 验证目录是否准确

3. **质量检查表**：
   - [ ] 所有 DOI 均通过 verify_itations.py 验证
   - [ ] 引文格式一致
   - [ ] 包括 PRISMA 流程图（用于系统评价）
   - [ ] 搜索方法完整记录
   - [ ] 明确规定纳入/排除标准
   - [ ] 结果按主题组织（不是逐项研究）
   - [ ] 质量评估已完成
   - [ ] 承认限制
   - [ ] 参考文献完整且准确
   - [ ] PDF 生成无错误

## 特定于数据库的搜索指南

### PubMed / PubMed 中心

通过`gget`技能访问：
```bash
# Search PubMed
gget search pubmed "CRISPR gene editing" -l 100

# Search with filters
# Use PubMed Advanced Search Builder to construct complex queries
# Then execute via gget or direct Entrez API
```

**搜索提示**：
- 使用 MeSH 术语：`"sickle cell disease"[MeSH]`
- 字段标签：`[Title]`、`[Title/Abstract]`、`[Author]`
- 日期过滤器：`2020:2024[Publication Date]`
- 布尔运算符：AND、OR、NOT
- 查看 MeSH 浏览器：https://meshb.nlm.nih.gov/search

### bioRxiv / medRxiv

通过`gget`技能访问：
```bash
gget search biorxiv "CRISPR sickle cell" -l 50
```

**重要注意事项**：
- 预印本未经同行评审
- 谨慎验证结果
- 检查预印本是否已发布（CrossRef）
- 注意预印本版本和日期

### arXiv

通过直接 API 或 WebFetch 访问：
```python
# Example search categories:
# q-bio.QM (Quantitative Methods)
# q-bio.GN (Genomics)
# q-bio.MN (Molecular Networks)
# cs.LG (Machine Learning)
# stat.ML (Machine Learning Statistics)

# Search format: category AND terms
search_query = "cat:q-bio.QM AND ti:\"single cell sequencing\""
```

### 语义学者

通过直接 API 访问（需要 API 密钥，或使用免费套餐）：
- 跨所有领域的 2 亿多篇论文
- 非常适合跨学科搜索
- 提供引文图表和论文推荐
- 用于查找极具影响力的论文

### 专业生物医学数据库

使用适当的技能：
- **ChEMBL**：`bioservices` 化学生物活性技能
- **UniProt**：`gget` 或 `bioservices` 蛋白质信息技能
- **KEGG**：`bioservices` 通路和基因技能
- **COSMIC**：`gget` 癌症突变技能
- **AlphaFold**：`gget alphafold` 用于蛋白质结构
- **PDB**：`gget` 或用于实验结构的直接 API

### 引文链接

通过引文网络扩展搜索：

1. **前向引用**（引用关键论文的论文）：
   - 使用谷歌学术搜索“被引用”
   - 使用 Semantic Scholar 或 OpenAlex API
   - 确定基于开创性工作的更新研究

2. **反向引用**（关键论文的参考文献）：
   - 从包含的论文中提取参考文献
   - 确定被高度引用的基础工作
   - 查找多项纳入研究引用的论文

## 引文风格指南

详细的格式指南位于`references/citation_styles.md`。快速参考：

### APA（第七版）
- 文内：（Smith 等人，2023）
- 参考文献：Smith, J. D.、Johnson, M. L. 和 Williams, K. R. (2023)。标题。 *期刊*，*22*(4)，301-318。 https://doi.org/10.xxx/yyy

### 自然
- 文本内：上标数字^1,2^
- 参考文献：Smith, J. D.、Johnson, M. L. 和 Williams, K. R. 标题。 *纳特。修订版药物发现* **22**, 301-318 (2023)。

### 温哥华
- 文本内：上标数字^1,2^
- 参考文献：Smith JD、Johnson ML、Williams KR。标题。 Nat Rev 药物发现。 2023；22(4):301-18。

**在最终确定之前，请务必使用 verify_itations.py 验证引文**。

## 最佳实践

### 搜索策略
1. **使用多个数据库**（至少3个）：确保全面覆盖
2. **包括预印本服务器**：捕获最新未发表的发现
3. **记录一切**：搜索字符串、日期、结果计数以确保可重复性
4. **测试和完善**：运行试点搜索、审查结果、调整搜索词

### 筛选和选择
1. **使用明确的标准**：筛选前记录纳入/排除标准
2. **系统筛选**：标题→摘要→全文
3. **文件排除**：记录排除研究的原因
4. **考虑双重筛选**：对于系统评价，请两名审稿人独立筛选

### 合成
1. **按主题组织**：按主题分组，而不是按个别研究分组
2. **跨研究综合**：比较、对比、识别模式
3. **保持批判性**：评估证据的质量和一致性
4. **找出差距**：记下遗漏或未充分研究的内容

### 质量和再现性
1. **评估研究质量**：使用适当的质量评估工具
2. **验证所有引文**：运行 verify_itations.py 脚本
3. **文档方法**：提供足够的细节供其他人复制
4. **遵循指南**：使用 PRISMA 进行系统评价

### 写作
1. **客观**：公平地提供证据，承认局限性
2. **系统化**：遵循结构化模板
3. **具体**：包括数字、统计数据、效应大小（如果有）
4. **清晰**：使用清晰的标题、逻辑流程、主题组织

## 要避免的常见陷阱

1. **单一数据库检索**：遗漏相关论文；总是搜索多个数据库
2. **没有搜索文档**：使得审查不可重复；记录所有搜索
3. **逐项研究总结**：缺乏综合；相反，按主题组织
4. **未经验证的引用**：导致错误；始终运行 verify_itations.py
5. **搜索范围太广**：产生数千个不相关的结果；用具体条款细化
6. **搜索范围太窄**：错过相关论文；包括同义词和相关术语
7. **忽略预印本**：错过最新发现；包括bioRxiv、medRxiv、arXiv
8. **无质量评估**：平等对待所有证据；评估和报告质量
9. **发表偏倚**：仅发表正面结果；注意潜在的偏见
10. **过时的搜索**：领域发展迅速；明确注明检索日期

## 工作流程示例

生物医学文献综述的完整工作流程：

```bash
# 1. Create review document from template
cp assets/review_template.md crispr_sickle_cell_review.md

# 2. Search multiple databases using appropriate skills
# - Use gget skill for PubMed, bioRxiv
# - Use direct API access for arXiv, Semantic Scholar
# - Export results in JSON format

# 3. Aggregate and process results
python scripts/search_databases.py combined_results.json \
  --deduplicate \
  --rank citations \
  --year-start 2015 \
  --year-end 2024 \
  --format markdown \
  --output search_results.md \
  --summary

# 4. Screen results and extract data
# - Manually screen titles, abstracts, full texts
# - Extract key data into the review document
# - Organize by themes

# 5. Write the review following template structure
# - Introduction with clear objectives
# - Detailed methodology section
# - Results organized thematically
# - Critical discussion
# - Clear conclusions

# 6. Verify all citations
python scripts/verify_citations.py crispr_sickle_cell_review.md

# Review the citation report
cat crispr_sickle_cell_review_citation_report.json

# Fix any failed citations and re-verify
python scripts/verify_citations.py crispr_sickle_cell_review.md

# 7. Generate professional PDF
python scripts/generate_pdf.py crispr_sickle_cell_review.md \
  --citation-style nature \
  --output crispr_sickle_cell_review.pdf

# 8. Review final PDF and markdown outputs
```

## 与其他技能的整合

这项技能与其他科学技能无缝配合：

### 数据库访问技巧
- **gget**：PubMed、bioRxiv、COSMIC、AlphaFold、Ensembl、UniProt
- **生物服务**：ChEMBL、KEGG、Reactome、UniProt、PubChem
- **datacommons-client**：人口统计、经济、健康统计

### 分析技巧
- **pydeseq2**：RNA-seq差异表达（用于方法部分）
- **scanpy**：单细胞分析（用于方法部分）
- **anndata**：单细胞数据（用于方法部分）
- **biopython**：序列分析（用于背景部分）

### 可视化技能
- **matplotlib**：生成图形和绘图以供审核
- **seaborn**：统计可视化

### 写作技巧
- **品牌指南**：将机构品牌应用于 PDF
- **内部通讯**：针对不同受众调整评论

## 资源

### 捆绑资源

**脚本：**
- `scripts/verify_citations.py`：验证 DOI 并生成格式化引文
- `scripts/generate_pdf.py`：将 Markdown 转换为专业 PDF
- `scripts/search_databases.py`：处理、删除重复数据并格式化搜索结果

**参考资料：**
- `references/citation_styles.md`：详细的引文格式指南（APA、Nature、温哥华、芝加哥、IEEE）
- `references/database_strategies.md`：全面的数据库搜索策略

**资产：**
- `assets/review_template.md`：包含所有部分的完整文献综述模板

### 外部资源

**指南：**
- PRISMA（系统评论）：http://www.prisma-statement.org/
- Cochrane 手册：https://training.cochrane.org/handbook
- AMSTAR 2（评论质量）：https://amstar.ca/

**工具：**
- MeSH 浏览器：https://meshb.nlm.nih.gov/search
- PubMed 高级搜索：https://pubmed.ncbi.nlm.nih.gov/advanced/
- 布尔搜索指南：https://www.ncbi.nlm.nih.gov/books/NBK3827/

**引文样式：**
- APA 风格：https://apastyle.apa.org/
- 自然作品集：https://www.nature.com/nature-portfolio/editorial-policies/reporting-standards
- NLM/温哥华：https://www.nlm.nih.gov/bsd/uniform_requirements.html

## 依赖关系

### 所需的 Python 包
```bash
uv pip install requests  # For citation verification
```

### 所需的系统工具
```bash
# For PDF generation
brew install pandoc  # macOS
apt-get install pandoc  # Linux

# For LaTeX (PDF generation)
brew install --cask mactex  # macOS
apt-get install texlive-xetex  # Linux
```

检查依赖关系：
```bash
python scripts/generate_pdf.py --check-deps
```

## 总结

这种文献综述技巧提供：

1. **系统方法**遵循学术最佳实践
2. **通过现有的科学技能进行多数据库集成**
3. **引文验证** 确保准确性和可信度
4. **专业输出** Markdown 和 PDF 格式
5.**全面指导**涵盖整个审核流程
6. **质量保证** 使用验证和确认工具
7. 通过详细的文档要求实现**可重复性**

进行全面、严格的文献综述，以满足学术标准，并提供任何领域当前知识的全面综合。