<!-- 此文件由机器翻译自 descriptors_viz.md -->

# Datamol 描述符和可视化参考

## 描述符模块 (`datamol.descriptors`)

描述符模块提供了用于计算分子属性和描述符的工具。

### 专用描述符函数

#### `dm.descriptors.n_aromatic_atoms(mol)`
计算芳香原子的数量。
- **返回**：整数计数
- **用例**：芳香度分析

#### `dm.descriptors.n_aromatic_atoms_proportion(mol)`
计算芳香原子与总重原子的比率。
- **返回**：在 0 和 1 之间浮动
- **用例**：量化芳香特征

#### `dm.descriptors.n_charged_atoms(mol)`
计算具有非零形式电荷的原子。
- **返回**：整数计数
- **用例**：电荷分布分析

#### `dm.descriptors.n_rigid_bonds(mol)`
计算不可旋转键（既不是单键也不是环键）。
- **返回**：整数计数
- **用例**：分子灵活性评估

#### `dm.descriptors.n_stereo_centers(mol)`
计算立体中心（手性中心）。
- **返回**：整数计数
- **用例**：立体化学分析

#### `dm.descriptors.n_stereo_centers_unspecified(mol)`
计算缺乏立体化学规范的立体中心。
- **返回**：整数计数
- **用例**：识别不完整的立体化学

### 批量描述符计算

#### `dm.descriptors.compute_many_descriptors(mol, properties_fn=None, add_properties=True)`
计算单个分子的多个分子特性。
- **参数**：
  - `properties_fn`：描述符函数的自定义列表
  - `add_properties`：包括其他计算属性
- **返回**：描述符名称→值对的字典
- **默认描述符包括**：
  - 分子量、LogP、氢键供体/受体数量
  - 芳香原子、立构中心、可旋转键
  - TPSA（拓扑极性表面积）
  - 环计数、杂原子计数
- **示例**：
  ```python
  mol = dm.to_mol("CCO")
  descriptors = dm.descriptors.compute_many_descriptors(mol)
  # Returns: {'mw': 46.07, 'logp': -0.03, 'hbd': 1, 'hba': 1, ...}
  ```

#### `dm.descriptors.batch_compute_many_descriptors(mols, properties_fn=None, add_properties=True, n_jobs=1, batch_size=None, progress=False)`
并行计算多个分子的描述符。
- **参数**：
  - `mols`：分子列表
  - `n_jobs`：并行作业数量（所有核心为-1）
  - `batch_size`：并行处理的块大小
  - `progress`：显示进度条
- **返回**：Pandas DataFrame，每个分子一行
- **示例**：
  <<<代码块_1>>>

### RDKit 描述符访问

#### `dm.descriptors.any_rdkit_descriptor(name)`
按名称从 RDKit 检索任何描述符函数。
- **参数**：`name` - 描述符函数名称（例如，“MolWt”、“TPSA”）
- **返回**：RDKit描述符函数
- **可用描述符**：来自 `rdkit.Chem.Descriptors` 和 `rdkit.Chem.rdMolDescriptors`
- **示例**：
  <<<代码块_2>>>

### 常见用例

**药物相似性过滤（Lipinski 的五法则）**：
<<<代码块_3>>>

**ADME性能分析**：
<<<代码块_4>>>

---

## 可视化模块 (`datamol.viz`)

viz 模块提供了将分子和构象异构体渲染为图像的工具。

### 主要可视化功能

#### `dm.viz.to_image(mols, legends=None, n_cols=4, use_svg=False, mol_size=(200, 200), highlight_atom=None, highlight_bond=None, outfile=None, max_mols=None, copy=True, indices=False, ...)`
从分子生成图像网格。
- **参数**：
  - `mols`：单个分子或分子列表
  - `legends`：作为标签的字符串或字符串列表（每个分子一个）
  - `n_cols`：每行的分子数（默认值：4）
  - `use_svg`：输出 SVG 格式（True）或 PNG（False，默认）
  - `mol_size`：方形图像的元组（宽度、高度）或单个整数
  - `highlight_atom`：要突出显示的原子索引（列表或字典）
  - `highlight_bond`：要突出显示的键索引（列表或字典）
  - `outfile`：保存路径（本地或远程，支持fsspec）
  - `max_mols`：要显示的最大分子数
  - `indices`：在结构上绘制原子索引（默认值：False）
  - `align`：使用 MCS（最大公共子结构）对齐分子
- **返回**：图像对象（可以在Jupyter中显示）或保存到文件
- **示例**：
  <<<代码块_5>>>

### 符合者可视化

#### `dm.viz.conformers(mol, n_confs=None, align_conf=True, n_cols=3, sync_views=True, remove_hs=True, ...)`
在网格布局中显示多个构象异构体。
- **参数**：
  - `mol`：具有嵌入构象异构体的分子
  - `n_confs`：要显示的构象异构体索引的数量或列表（无 = 全部）
  - `align_conf`：对齐构象异构体进行比较（默认值：True）
  - `n_cols`：网格列（默认值：3）
  - `sync_views`：交互时同步 3D 视图（默认值：True）
  - `remove_hs`：为了清晰起见，删除氢（默认值：True）
- **返回**：符合异构体可视化网格
- **用例**：比较构象多样性
- **示例**：
  <<<代码块_6>>>

### 圆形网格可视化
#### `dm.viz.circle_grid(center_mol, circle_mols, mol_size=200, circle_margin=50, act_mapper=None, ...)`
使用中心分子创建同心环可视化。
- **参数**：
  - `center_mol`：中心分子
  - `circle_mols`：分子列表列表（每个环一个列表）
  - `mol_size`：每个分子的图像大小
  - `circle_margin`：环之间的间距（默认值：50）
  - `act_mapper`：用于颜色编码的活动映射字典
- **返回**：圆形网格图像
- **用例**：可视化分子邻域、SAR 分析、相似性网络
- **示例**：
  ```python
  # Show a reference molecule surrounded by similar compounds
  dm.viz.circle_grid(
      center_mol=reference,
      circle_mols=[nearest_neighbors, second_tier]
  )
  ```

### 可视化最佳实践

1. **为了清晰起见，使用图例**：始终使用 SMILES、ID 或活性值来标记分子
2. **对齐相关分子**：使用`to_image()`中的`align=True`进行SAR分析
3. **调整网格大小**：根据分子数和显示宽度设置`n_cols`
4. **使用 SVG 进行出版物**：为可缩放矢量图形设置 `use_svg=True`
5. **突出子结构**：使用`highlight_atom`和`highlight_bond`来强调功能
6. **保存大网格**：使用`outfile`参数保存而不是在内存中显示