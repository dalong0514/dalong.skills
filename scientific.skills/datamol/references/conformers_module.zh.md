<!-- 此文件由机器翻译自 conformers_module.md -->

# Datamol Conformers 模块参考

`datamol.conformers` 模块提供用于生成和分析 3D 分子构象的工具。

## 顺应者一代

### `dm.conformers.generate(mol, n_confs=None, rms_cutoff=None, minimize_energy=True, method='ETKDGv3', add_hs=True, ...)`
生成 3D 分子构象异构体。
- **参数**：
  - `mol`：输入分子
  - `n_confs`：要生成的构象异构体数量（如果没有，则根据可旋转键自动确定）
  - `rms_cutoff`：以 Ångströms 为单位的 RMS 阈值，用于过滤相似的构象异构体（删除重复项）
  - `minimize_energy`：应用 UFF 能量最小化（默认值：True）
  - `method`：嵌入方法 - 选项：
    - `'ETDG'` - 实验扭转距离几何
    - `'ETKDG'` - ETDG 以及额外的基础知识
    - `'ETKDGv2'` - 增强版本 2
    - `'ETKDGv3'` - 增强版本 3（默认，推荐）
  - `add_hs`：嵌入前添加氢（默认值：True，对质量至关重要）
  - `random_seed`：设置可重复性
- **返回**：具有嵌入构象异构体的分子
- **示例**：
  ```python
  mol = dm.to_mol("CCO")
  mol_3d = dm.conformers.generate(mol, n_confs=10, rms_cutoff=0.5)
  conformers = mol_3d.GetConformers()  # Access all conformers
  ```

## 一致性聚类

### `dm.conformers.cluster(mol, rms_cutoff=1.0, already_aligned=False, centroids=False)`
按 RMS 距离对构象异构体进行分组。
- **参数**：
  - `rms_cutoff`：以 Ångströms 为单位的聚类阈值（默认值：1.0）
  - `already_aligned`：构象异构体是否预对齐
  - `centroids`：返回质心构象异构体（True）或簇组（False）
- **返回**：簇信息或质心构象异构体
- **用例**：识别不同的构象家族

### `dm.conformers.return_centroids(mol, conf_clusters, centroids=True)`
从簇中提取代表性构象异构体。
- **参数**：
  - `conf_clusters`：来自 `cluster()` 的簇索引序列
  - `centroids`：返回单个分子（True）或分子列表（False）
- **返回**：质心构象异构体

## 一致性分析

### `dm.conformers.rmsd(mol)`
计算所有构象异构体的成对 RMSD 矩阵。
- **要求**：至少 2 个符合者
- **返回**：RMSD 值的 NxN 矩阵
- **用例**：量化构象多样性

### `dm.conformers.sasa(mol, n_jobs=1, ...)`
使用 FreeSASA 计算溶剂可及表面积 (SASA)。
- **参数**：
  - `n_jobs`：多个构象异构体的并行化
- **返回**：SASA 值数组（每个符合者一个）
- **存储**：作为属性 `'rdkit_free_sasa'` 存储在每个构象异构体中的值
- **示例**：
  <<<代码块_1>>>

## 低级一致性操作

### `dm.conformers.center_of_mass(mol, conf_id=-1, use_atoms=True, round_coord=None)`
计算分子中心。
- **参数**：
  - `conf_id`：符合者索引（第一个符合者为-1）
  - `use_atoms`：使用原子质量（True）或几何中心（False）
  - `round_coord`：舍入的小数精度
- **返回**：中心的 3D 坐标
- **用例**：将分子居中以进行可视化或对齐

### `dm.conformers.get_coords(mol, conf_id=-1)`
从构象异构体中检索原子坐标。
- **返回**：原子位置的 Nx3 numpy 数组
- **示例**：
  <<<代码块_2>>>

### `dm.conformers.translate(mol, conf_id=-1, transform_matrix=None)`
使用变换矩阵重新定位构象异构体。
- **修改**：就地操作
- **用例**：对齐或重新定位分子

## 工作流程示例

<<<代码块_3>>>

## 关键概念

- **距离几何**：根据连通性信息生成 3D 结构的方法
- **ETKDG**：使用实验扭转角偏好和附加化学知识
- **RMS Cutoff**：较低的值 = 更多独特的构象异构体；更高的值=更少、更独特的构象异构体
- **能量最小化**：将结构放松到最接近的局部能量最小值
- **氢**：对于精确的 3D 几何形状至关重要 - 在嵌入过程中始终包含