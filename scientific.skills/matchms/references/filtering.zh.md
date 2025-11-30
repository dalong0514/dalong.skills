<!-- 此文件由机器翻译自 filtering.md -->

# Matchms 过滤函数参考

本文档提供了 matchms 中用于处理质谱数据的所有过滤功能的综合参考。

## 元数据处理过滤器

### 化合物和化学信息

**添加化合物名称（光谱）**
- 将化合物名称添加到正确的元数据字段
- 标准化化合物名称存储位置

**clean_compound_name（光谱）**
- 从化合物名称中删除常见的不需要的添加内容
- 清理格式不一致的问题

**derive_adduct_from_name（谱）**
- 从化合物名称中提取加合物信息
- 将加合物符号移至正确的元数据字段

**从名称导出公式（频谱）**
- 检测化合物名称中的化学式
- 将公式重新定位到适当的元数据字段

**从化合物名称派生注释（光谱）**
- 使用化合物名称从 PubChem 检索 SMILES/InChI
- 自动注释化学结构

### 化学结构转换

**derive_inchi_from_smiles（光谱）**
- 从 SMILES 字符串生成 InChI
- 需要rdkit库

**derive_inchikey_from_inchi（频谱）**
- 从 InChI 计算 InChIKey
- 27 个字符的哈希标识符

**derive_smiles_from_inchi（光谱）**
- 从 InChI 表示创建 SMILES
- 需要rdkit库

**repair_inchi_inchikey_smiles（光谱）**
- 更正错误放置的化学标识符
- 修复元数据字段混乱

**repair_not_matching_annotation（光谱）**
- 确保 SMILES、InChI 和 InChIKey 之间的一致性
- 验证化学结构注释匹配

**add_fingerprint（光谱，fingerprint_type =“日光”，nbits = 2048，半径= 2）**
- 生成分子指纹以进行相似性计算
- 指纹类型：“日光”、“摩根1”、“摩根2”、“摩根3”
- 与指纹相似度评分一起使用

### 质量和电荷信息

**add_precursor_mz（谱）**
- 标准化前体 m/z 值
- 标准化前体质量元数据

**add_parent_mass（光谱，estimate_from_adduct=True）**
- 根据前体 m/z 和加合物计算中性母体质量
- 如果不能直接获得，可以根据加合物进行估计

**正确电荷（频谱）**
- 将电荷值与离子模式对齐
- 确保电荷符号与电离模式匹配

**make_charge_int（频谱）**
- 将电荷转换为整数格式
- 标准化电荷表示

**clean_adduct（光谱）**
- 标准化加合物符号
- 纠正常见的加合物格式问题

**解释_pepmass（光谱）**
- 将 pepmass 字段解析为组件值
- 从组合场中提取母体 m/z 和强度

### 离子模式和验证

**derive_ionmode（频谱）**
- 根据加合物信息确定离子模式
- 从加合物类型推断正/负模式

**require_ Correct_ionmode（光谱，ion_mode）**
- 按指定离子模式过滤光谱
- 如果 ionmode 不匹配则返回 None
- 使用：`spectrum = require_correct_ionmode(spectrum, "positive")`

**require_precursor_mz(频谱，minimum_accepted_mz=0.0)**
- 验证前体 m/z 的存在和值
- 如果缺失或低于阈值则返回 None

**require_precursor_below_mz(频谱，maximum_accepted_mz=1000.0)**
- 强制执行最大母离子 m/z 限制
- 如果前体超过阈值则返回 None

### 保留信息

**添加保留时间（光谱）**
- 将保留时间协调为浮点值
- 标准化 RT 元数据字段

**添加保留指数（谱）**
- 在标准化字段中存储保留指数
- 规范化 RI 元数据

### 数据协调

**harmonize_undefine_inchi（光谱，未定义=“”，别名=无）**
- 标准化未定义/空的 InChI 条目
- 用一致的值替换各种“未知”的表示

**harmonize_undefined_inchikey（光谱，未定义=“”，别名=无）**
- 标准化未定义/空的 InChIKey 条目
- 统一缺失数据表示

**harmonize_undefined_smiles（光谱，未定义=“”，别名=无）**
- 标准化未定义/空的 SMILES 条目
- 一致处理缺失的结构数据

### 维修和质量职能

**repair_adduct_based_on_smiles（光谱，mass_tolerance=0.1）**
- 使用 SMILES 和质量匹配校正加合物
- 验证加合物与计算质量的匹配

**repair_parent_mass_is_mol_wt（光谱，mass_tolerance=0.1）**
- 将分子量转换为单同位素质量
- 修复常见的元数据混乱

**repair_precursor_is_parent_mass（光谱）**
- 修复交换的前体/母体质量值
- 纠正字段错误分配

**repair_smiles_of_salts（光谱，mass_tolerance=0.1）**
- 去除盐成分以匹配母体质量
- 提取相关分子片段

**require_parent_mass_match_smiles（光谱，mass_tolerance=0.1）**
- 根据 SMILES 计算的质量验证母体质量
- 如果质量在公差范围内不匹配，则返回“无”

**require_valid_annotation（频谱）**
- 确保化学注释完整、一致
- 验证 SMILES、InChI 和 InChIKey 的存在和一致性

## 峰值处理滤波器

### 标准化和选择

**归一化强度（频谱）**
- 将峰值强度缩放至单位高度（最大值 = 1.0）
- 相似性计算的基本预处理步骤

**按强度选择（光谱，intensity_from=0.0，intensity_to=1.0）**
- 保留指定绝对强度范围内的峰值
- 按原始强度值过滤

**按相对强度选择（光谱，intensity_from=0.0，intensity_to=1.0）**
- 将峰保持在相对强度范围内
- 过滤器作为最大强度的一部分

**select_by_mz（光谱，mz_from=0.0，mz_to=1000.0）**
- 按 m/z 值范围过滤峰值
- 删除指定 m/z 窗口之外的峰

### 峰值降低和过滤

**reduce_to_number_of_peaks（光谱，n_max=无，ratio_desired=无）**
- 超过最大值时删除最低强度峰值
- 可以指定绝对数量或比例
- 使用：`spectrum = reduce_to_number_of_peaks(spectrum, n_max=100)`

**remove_peaks_around_precursor_mz（光谱，mz_tolerance=17）**
- 消除母体耐受范围内的峰
- 去除前体和同位素峰
- 基于片段的相似性的通用预处理

**remove_peaks_outside_top_k（频谱，k=10，ratio_desired=None）**
- 仅保留 k 个最高强度峰附近的峰
- 专注于最具信息性的信号

**require_minimum_number_of_peaks（光谱，n_required=10）**
- 丢弃峰值不足的光谱
- 质量控制过滤器
- 如果峰值计数低于阈值则返回 None

**require_minimum_number_of_high_peaks（光谱，n_required=5，intensity_threshold=0.05）**
- 删除缺乏高强度峰的光谱
- 确保数据质量
- 如果没有足够的峰值高于阈值，则返回 None

### 损失计算

**add_losses（频谱，loss_mz_from=5.0，loss_mz_to=200.0）**
- 从前体质量中得出中性损失
- 计算损失=前体_mz - 片段_mz
- 为 NeutralLossesCosine 评分添加损失谱

## 管道函数

**默认过滤器（频谱）**
- 按顺序应用九个基本元数据过滤器：
  1.make_charge_int
  2.add_precursor_mz
  3.添加保留时间
  4.添加保留索引
  5. 从名称中派生加合物
  6. 从名称导出公式
  7. 清洁化合物名称
  8. 和谐未定义微笑
  9.协调_未定义_英寸
- 元数据协调的推荐起点

**频谱处理器（滤波器）**
- 协调多过滤器管道
- 接受过滤功能列表
- 示例：
```python
from matchms import SpectrumProcessor
processor = SpectrumProcessor([
    default_filters,
    normalize_intensities,
    lambda s: select_by_relative_intensity(s, intensity_from=0.01)
])
processed = processor(spectrum)
```

## 常见过滤器组合

### 标准预处理管道
<<<代码块_1>>>

### 质量控制流程
<<<代码块_2>>>

### 化学注释管道
<<<代码块_3>>>

### 高峰清洁管道
<<<代码块_4>>>

## 过滤器使用注意事项

1. **顺序很重要**：按逻辑顺序应用过滤器（例如，在相对强度选择之前进行归一化）
2. **过滤器返回 None**：许多过滤器对于无效光谱返回 None；在继续之前检查 None
3. **不变性**：过滤器通常返回修改后的副本；将结果重新分配给变量
4. **管道效率**：使用SpectrumProcessor进行一致的多频谱处理
5. **文档**：详细参数请参见matchms.readthedocs.io/en/latest/api/matchms.filtering.html