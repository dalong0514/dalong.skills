<!-- 此文件由机器翻译自 protein_optimization.md -->

# 蛋白质序列优化

## 概述

在提交蛋白质序列进行实验测试之前，使用计算工具优化序列以提高表达、溶解度和稳定性。这种预筛选降低了实验成本并提高了成功率。

## 常见的蛋白质表达问题

### 1. 未配对的半胱氨酸

**问题：**
- 未配对的半胱氨酸形成不需要的二硫键
- 导致聚集和错误折叠
- 降低表达产量和稳定性

**解决方案：**
- 除非功能需要，否则删除未配对的半胱氨酸
- 将半胱氨酸适当配对为结构二硫键
- 在非关键位置替换为丝氨酸或丙氨酸

**示例：**
```python
# Check for cysteine pairs
from Bio.Seq import Seq

def check_cysteines(sequence):
    cys_count = sequence.count('C')
    if cys_count % 2 != 0:
        print(f"Warning: Odd number of cysteines ({cys_count})")
    return cys_count
```

### 2. 过度疏水性

**问题：**
- 长疏水斑块促进聚集
- 暴露的疏水残基导致蛋白质聚集
- 在水性缓冲液中溶解度差

**解决方案：**
- 保持平衡的水疗特性
- 在域之间使用短而灵活的链接器
- 减少表面暴露的疏水性残留物

**指标：**
- Kyte-Doolittle 水病图
- GRAVY 评分（水病总平均值）
- pSAE（溶剂可及的疏水残基百分比）

### 3. 低溶解度

**问题：**
- 表达或纯化过程中蛋白质沉淀
- 包涵体形成
- 下游加工困难

**解决方案：**
- 使用溶解度预测工具进行预筛选
- 应用序列优化算法
- 如果需要，添加增溶标签

## 用于优化的计算工具

### NetSolP - 初始溶解度筛选

**目的：** 过滤序列的快速溶解度预测。

**方法：** 在大肠杆菌表达数据上训练的机器学习模型。

**用途：**
<<<代码块_1>>>

**释义：**
- 分数 > 0.5：可能溶解
- 分数 < 0.5：可能不溶
- 用于在更昂贵的预测之前进行初始过滤

**何时使用：**
- 大型库的首次过滤
- 快速验证设计的序列
- 确定实验测试序列的优先级

### SoluProt - 综合溶解度预测

**目的：** 具有更高准确度的高级溶解度预测。

**方法：** 结合序列和结构特征的深度学习模型。

**用途：**
<<<代码块_2>>>

**释义：**
- 分数 > 0.6：高溶解度置信度
- 分数 0.4-0.6：不确定，可能需要优化
- 分数 < 0.4：可能有问题

**何时使用：**
- 初始 NetSolP 过滤后
- 当需要更高的预测精度时
- 在进行昂贵的合成/测试之前

### SolubleMPNN - 序列重新设计

**目的：** 重新设计蛋白质序列以提高溶解度，同时保持功能。

**方法：** 图神经网络建议突变以增加溶解度。

**用途：**
<<<代码块_3>>>

**设计策略：**
- **保守**（温度=0.1）：变化最小，更安全
- **中等**（温度=0.3）：变化与安全之间的平衡
- **攻击性**（温度=0.5）：更多突变，更高风险

**何时使用：**
- 序列优化的主要工具
- 改进有问题的序列的默认起点
- 生成多种可溶性变体

**最佳实践：**
- 每个序列生成 10-50 个变体
- 使用可用的结构信息（提高准确性）
- 验证关键功能残基是否得到保留
- 测试多种温度设置

### ESM（进化尺度建模）- 序列似然

**目的：** 根据进化模式评估蛋白质序列的“自然”程度。

**方法：** 在数百万个自然序列上训练的蛋白质语言模型。

**用途：**
<<<代码块_4>>>

**释义：**
- 更高的分数→更“自然”的序列
- 用于避免不太可能发生的突变
- 与功能需求的平衡

**何时使用：**
- 过滤合成设计
- 比较 SolubleMPNN 变体
- 确保序列不会太人为
- 避免表达瓶颈

**与设计整合：**
<<<代码块_5>>>

### ipTM - 界面稳定性（AlphaFold-Multimer）

**目的：** 评估蛋白质-蛋白质界面稳定性和结合置信度。

**方法：** 界面根据 AlphaFold-Multimer 预测预测 TM 分数。

**用途：**
<<<代码块_6>>>

**释义：**
- ipTM > 0.7：强预测接口
- ipTM 0.5-0.7：中等接口置信度
- ipTM < 0.5：接口弱，考虑重新设计

**何时使用：**
- 抗体-抗原设计
- 蛋白质-蛋白质相互作用工程
- 验证绑定接口
- 比较界面变体

### pSAE - 溶剂可及的疏水性残留物

**目的：** 量化促进聚集的暴露疏水残基。

**方法：** 计算疏水残基占据的溶剂可及表面积 (SASA) 的百分比。

**用途：**
```python
# Requires structure (PDB file or AlphaFold prediction)
# Install: uv pip install biopython

from Bio.PDB import PDBParser, DSSP
import numpy as np

def calculate_psae(pdb_file):
    """
    Calculate percent Solvent-Accessible hydrophobic residues (pSAE)

    Lower pSAE = better solubility
    """

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    # Run DSSP to get solvent accessibility
    model = structure[0]
    dssp = DSSP(model, pdb_file, acc_array='Wilke')

    hydrophobic = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO']

    total_sasa = 0
    hydrophobic_sasa = 0

    for residue in dssp:
        res_name = residue[1]
        rel_accessibility = residue[3]

        total_sasa += rel_accessibility
        if res_name in hydrophobic:
            hydrophobic_sasa += rel_accessibility

    psae = (hydrophobic_sasa / total_sasa) * 100

    return psae

# Example
pdb_file = "protein_structure.pdb"
psae_score = calculate_psae(pdb_file)
print(f"pSAE: {psae_score:.2f}%")

# Interpretation
if psae_score < 25:
    print("Good solubility expected")
elif psae_score < 35:
    print("Moderate solubility")
else:
    print("High aggregation risk")
```

**释义：**
- pSAE < 25%：聚集风险低
- pSAE 25-35%：中等风险
- pSAE > 35%：高聚集风险

**何时使用：**
- 分析设计结构
- AlphaFold 后验证
- 识别聚合热点
- 引导表面突变

## 推荐的优化工作流程

### 第 1 步：初步筛选（快速）

```python
def initial_screening(sequences):
    """
    Quick first-pass filtering using NetSolP
    Filters out obviously problematic sequences
    """
    passed = []
    for name, seq in sequences.items():
        netsolp_score = predict_solubility_netsolp(seq)
        if netsolp_score > 0.5:
            passed.append((name, seq))

    return passed
```

### 第 2 步：详细评估（中等）

```python
def detailed_assessment(filtered_sequences):
    """
    More thorough analysis with SoluProt and ESM
    Ranks sequences by multiple criteria
    """
    results = []
    for name, seq in filtered_sequences:
        soluprot_score = predict_solubility(seq)
        esm_score = score_sequence_esm(seq)

        combined_score = soluprot_score * 0.7 + esm_score * 0.3

        results.append({
            'name': name,
            'sequence': seq,
            'soluprot': soluprot_score,
            'esm': esm_score,
            'combined': combined_score
        })

    results.sort(key=lambda x: x['combined'], reverse=True)
    return results
```

### 步骤 3：序列优化（如果需要）

```python
def optimize_problematic_sequences(sequences_needing_optimization):
    """
    Use SolubleMPNN to redesign problematic sequences
    Returns improved variants
    """
    optimized = []
    for name, seq in sequences_needing_optimization:
        # Generate multiple variants
        variants = optimize_sequence(
            sequence=seq,
            num_variants=10,
            temperature=0.2
        )

        # Score variants with ESM
        for variant in variants:
            variant['esm_score'] = score_sequence_esm(variant['sequence'])

        # Keep best variants
        variants.sort(
            key=lambda x: x['solubility_score'] * x['esm_score'],
            reverse=True
        )

        optimized.extend(variants[:3])  # Top 3 variants per sequence

    return optimized
```

### 步骤 4：基于结构的验证（对于关键序列）

```python
def structure_validation(top_candidates):
    """
    Predict structures and calculate pSAE for top candidates
    Final validation before experimental testing
    """
    validated = []
    for candidate in top_candidates:
        # Predict structure with AlphaFold
        structure_pdb = predict_structure_alphafold(candidate['sequence'])

        # Calculate pSAE
        psae = calculate_psae(structure_pdb)

        candidate['psae'] = psae
        candidate['pass_structure_check'] = psae < 30

        validated.append(candidate)

    return validated
```

### 完整的工作流程示例

```python
def complete_optimization_pipeline(initial_sequences):
    """
    End-to-end optimization pipeline

    Input: Dictionary of {name: sequence}
    Output: Ranked list of optimized, validated sequences
    """

    print("Step 1: Initial screening with NetSolP...")
    filtered = initial_screening(initial_sequences)
    print(f"  Passed: {len(filtered)}/{len(initial_sequences)}")

    print("Step 2: Detailed assessment with SoluProt and ESM...")
    assessed = detailed_assessment(filtered)

    # Split into good and needs-optimization
    good_sequences = [s for s in assessed if s['soluprot'] > 0.6]
    needs_optimization = [s for s in assessed if s['soluprot'] <= 0.6]

    print(f"  Good sequences: {len(good_sequences)}")
    print(f"  Need optimization: {len(needs_optimization)}")

    if needs_optimization:
        print("Step 3: Optimizing problematic sequences with SolubleMPNN...")
        optimized = optimize_problematic_sequences(needs_optimization)
        all_sequences = good_sequences + optimized
    else:
        all_sequences = good_sequences

    print("Step 4: Structure-based validation for top candidates...")
    top_20 = all_sequences[:20]
    final_validated = structure_validation(top_20)

    # Final ranking
    final_validated.sort(
        key=lambda x: (
            x['pass_structure_check'],
            x['combined'],
            -x['psae']
        ),
        reverse=True
    )

    return final_validated

# Usage
initial_library = {
    'variant_1': 'MKVLWAALLGLLGAAA...',
    'variant_2': 'MATGVLWAALLGLLGA...',
    # ... more sequences
}

optimized_library = complete_optimization_pipeline(initial_library)

# Submit top sequences to Adaptyv
top_sequences_for_testing = optimized_library[:50]
```

## 最佳实践总结

1. **在实验测试之前始终进行预筛选**
2. **首先使用NetSolP**快速过滤大型库
3. **应用SolubleMPNN**作为默认优化工具
4. **使用 ESM 验证**以避免不自然的序列
5. **计算 pSAE** 以进行基于结构的验证
6. **每个设计测试多个变体**以考虑预测的不确定性
7. **保留对照** - 包括野生型或已知良好的序列
8. **迭代** - 使用实验结果来完善预测

## 与 Adaptyv 集成

计算优化后，将序列提交给 Adaptyv：

```python
# After optimization pipeline
optimized_sequences = complete_optimization_pipeline(initial_library)

# Prepare FASTA format
fasta_content = ""
for seq_data in optimized_sequences[:50]:  # Top 50
    fasta_content += f">{seq_data['name']}\n{seq_data['sequence']}\n"

# Submit to Adaptyv
import requests
response = requests.post(
    "https://kq5jp7qj7wdqklhsxmovkzn4l40obksv.lambda-url.eu-central-1.on.aws/experiments",
    headers={"Authorization": f"Bearer {api_key}"},
    json={
        "sequences": fasta_content,
        "experiment_type": "expression",
        "metadata": {
            "optimization_method": "SolubleMPNN_ESM_pipeline",
            "computational_scores": [s['combined'] for s in optimized_sequences[:50]]
        }
    }
)
```

## 故障排除

**问题：所有序列在溶解度预测方面得分都很差**
- 检查序列是否含有异常氨基酸
- 验证 FASTA 格式是否正确
- 考虑蛋白质家族是否天然具有低溶解度
- 尽管有预测，但可能需要实验验证

**问题：SolubleMPNN 改变了功能上重要的残基**
- 提供结构文件以保留空间约束
- 掩盖突变的关键残基
- 较低的温度参数以进行保守的改变
- 手动恢复有问题的突变

**问题：优化后ESM分数较低**
- 优化可能过于激进
- 在 SolubleMPNN 中尝试较低的温度
- 溶解性与自然性之间的平衡
- 考虑到某些优化可能需要非自然突变

**问题：预测与实验结果不符**
- 预测是概率性的，而不是确定性的
- 宿主系统和条件影响表达
- 某些蛋白质可能需要实验验证
- 使用预测作为丰富，而不是绝对过滤器