<!-- 此文件由机器翻译自 SKILL.md -->

---
名称：科学可视化
描述：“使用 matplotlib/seaborn/plotly 创建出版图。多面板布局、误差线、重要性标记、色盲安全、导出 PDF/EPS/TIFF，用于期刊就绪的科学绘图。”
---

# 科学可视化

## 概述

科学可视化将数据转化为清晰、准确的数字以供发布。使用多面板布局、误差线、重要性标记和色盲安全调色板创建日志就绪绘图。使用 matplotlib、seaborn 和plotly 将手稿导出为 PDF/EPS/TIFF。

## 何时使用此技能

该技能应该在以下情况下使用：
- 为科学手稿创建绘图或可视化
- 准备提交期刊的图表（《Nature》、《Science》、《Cell》、《PLOS》等）
- 确保人物对色盲者友好且易于理解
- 制作具有一致样式的多面板人物
- 以正确的分辨率和格式导出图形
- 遵循特定的出版指南
- 改进现有数据以满足出版标准
- 创建需要在彩色和灰度下工作的图形

## 快速入门指南

### 基本出版质量图

```python
import matplotlib.pyplot as plt
import numpy as np

# Apply publication style (from scripts/style_presets.py)
from style_presets import apply_publication_style
apply_publication_style('default')

# Create figure with appropriate size (single column = 3.5 inches)
fig, ax = plt.subplots(figsize=(3.5, 2.5))

# Plot data
x = np.linspace(0, 10, 100)
ax.plot(x, np.sin(x), label='sin(x)')
ax.plot(x, np.cos(x), label='cos(x)')

# Proper labeling with units
ax.set_xlabel('Time (seconds)')
ax.set_ylabel('Amplitude (mV)')
ax.legend(frameon=False)

# Remove unnecessary spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Save in publication formats (from scripts/figure_export.py)
from figure_export import save_publication_figure
save_publication_figure(fig, 'figure1', formats=['pdf', 'png'], dpi=300)
```

### 使用预配置样式

使用 `assets/` 中的 matplotlib 样式文件应用特定于期刊的样式：

<<<代码块_1>>>

### Seaborn 快速入门

对于统计图，请使用带有发布样式的seaborn：

<<<代码块_2>>>

## 核心原则和最佳实践

### 1. 分辨率和文件格式

**关键要求**（详细信息参见`references/publication_guidelines.md`）：
- **光栅图像**（照片、显微镜）：300-600 DPI
- **线条艺术**（图表、绘图）：600-1200 DPI 或矢量格式
- **矢量格式**（首选）：PDF、EPS、SVG
- **光栅格式**：TIFF、PNG（科学数据切勿使用 JPEG）

**实施：**
<<<代码块_3>>>

### 2. 颜色选择 - 色盲辅助功能

**始终使用色盲友好的调色板**（详细信息请参见`references/color_palettes.md`）：

**推荐：Okabe-Ito 调色板**（可区分所有类型的色盲）：
<<<代码块_4>>>

**对于热图/连续数据：**
- 使用感知统一的颜色图：`viridis`、`plasma`、`cividis`
- 避免红绿发散地图（使用 `PuOr`、`RdBu`、`BrBG` 代替）
- 切勿使用 `jet` 或 `rainbow` 颜色图

**始终以灰度测试图形**以确保可解释性。

### 3. 版式和文本

**字体指南**（详细信息参见`references/publication_guidelines.md`）：
- 无衬线字体：Arial、Helvetica、Calibri
- **最终打印尺寸**的最小尺寸：
  - 轴标签：7-9 pt
  - 刻度标签：6-8 pt
  - 面板标签：8-12 pt（粗体）
- 标签的句子大小写：“时间（小时）”而不是“时间（小时）”
- 始终在括号中包含单位

**实施：**
<<<代码块_5>>>

### 4. 人物尺寸

**特定于期刊的宽度**（详细信息请参见`references/journal_requirements.md`）：
- **自然**：单89毫米，双183毫米
- **科学**：单 55 毫米，双 175 毫米
- **电池**：单 85 毫米，双 178 毫米

**检查图形尺寸合规性：**
<<<代码块_6>>>

### 5. 多面板人物

**最佳实践：**
- 用粗体字母标记面板：**A**、**B**、**C**（大多数期刊为大写字母，Nature 为小写字母）
- 在所有面板上保持一致的样式
- 尽可能沿边缘对齐面板
- 在面板之间使用足够的空白

**示例实现**（完整代码请参见`references/matplotlib_examples.md`）：
```python
from string import ascii_uppercase

fig = plt.figure(figsize=(7, 4))
gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.4)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
# ... create other panels ...

# Add panel labels
for i, ax in enumerate([ax1, ax2, ...]):
    ax.text(-0.15, 1.05, ascii_uppercase[i], transform=ax.transAxes,
            fontsize=10, fontweight='bold', va='top')
```

## 常见任务

### 任务 1：创建可供发布的线图

有关完整代码，请参阅 `references/matplotlib_examples.md` 示例 1。

**关键步骤：**
1.应用出版风格
2.为目标期刊设置合适的图形尺寸
3.使用色盲友好的颜色
4. 添加具有正确表示的误差线（SEM、SD 或 CI）
5. 用单位标记轴
6.去除不必要的刺
7.保存为矢量格式

**使用seaborn进行自动置信区间：**
```python
import seaborn as sns
fig, ax = plt.subplots(figsize=(5, 3))
sns.lineplot(data=timeseries, x='time', y='measurement',
             hue='treatment', errorbar=('ci', 95), 
             markers=True, ax=ax)
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Measurement (AU)')
sns.despine()
```

### 任务 2：创建多面板图

有关完整代码，请参阅 `references/matplotlib_examples.md` 示例 2。

**关键步骤：**
1.使用`GridSpec`进行灵活布局
2. 确保跨面板的样式一致
3.添加粗体面板标签（A、B、C等）
4. 对齐相关面板
5. 验证所有文本在最终尺寸下均可读

### 任务 3：创建具有正确颜色图的热图

有关完整代码，请参阅 `references/matplotlib_examples.md` 示例 4。

**关键步骤：**
1. 使用感知统一的颜色图（`viridis`、`plasma`、`cividis`）
2. 包含带标签的颜色条
3. 对于发散数据，使用色盲安全发散图 (`RdBu_r`, `PuOr`)
4.为发散地图设置适当的中心值
5. 灰度外观测试

**使用seaborn作为相关矩阵：**
```python
import seaborn as sns
fig, ax = plt.subplots(figsize=(5, 4))
corr = df.corr()
mask = np.triu(np.ones_like(corr, dtype=bool))
sns.heatmap(corr, mask=mask, annot=True, fmt='.2f',
            cmap='RdBu_r', center=0, square=True,
            linewidths=1, cbar_kws={'shrink': 0.8}, ax=ax)
```

### 任务 4：为特定期刊准备图表

**工作流程：**
1. 检查期刊要求：`references/journal_requirements.md`
2.为journal配置matplotlib：
   ```python
   from style_presets import configure_for_journal
   configure_for_journal('nature', figure_width='single')
   ```
3. 创建图形（将自动正确调整大小）
4. 导出日志规格：
   ```python
   from figure_export import save_for_journal
   save_for_journal(fig, 'figure1', journal='nature', figure_type='line_art')
   ```

### 任务 5：修复现有图形以满足出版标准

**清单方法**（完整清单位于`references/publication_guidelines.md`中）：

1. **检查分辨率**：验证 DPI 是否满足期刊要求
2. **检查文件格式**：绘图使用矢量，图像使用 TIFF/PNG
3. **检查颜色**：确保色盲友好
4. **检查字体**：最终尺寸最小 6-7 pt，无衬线
5. **检查标签**：所有标有单位的轴
6. **检查尺寸**：匹配日记帐栏宽度
7. **测试灰度**：无需颜色即可解释的图
8. **删除图表垃圾**：没有不必要的网格、3D效果、阴影

### 任务 6：创建色盲友好的可视化效果

**策略：**
1. 使用 `assets/color_palettes.py` 中认可的调色板
2.添加冗余编码（线条样式、标记、图案）
3.用色盲模拟器测试
4、保证灰度兼容

**示例：**
```python
from color_palettes import apply_palette
import matplotlib.pyplot as plt

apply_palette('okabe_ito')

# Add redundant encoding beyond color
line_styles = ['-', '--', '-.', ':']
markers = ['o', 's', '^', 'v']

for i, (data, label) in enumerate(datasets):
    plt.plot(x, data, linestyle=line_styles[i % 4],
             marker=markers[i % 4], label=label)
```

## 统计严谨性

**始终包括：**
- 误差线（SD、SEM 或 CI - 在标题中指定）
- 图表或标题中的样本大小（n）
- 统计显着性标记（*、**、***）
- 尽可能单独的数据点（不仅仅是汇总统计数据）

**统计示例：**
```python
# Show individual points with summary statistics
ax.scatter(x_jittered, individual_points, alpha=0.4, s=8)
ax.errorbar(x, means, yerr=sems, fmt='o', capsize=3)

# Mark significance
ax.text(1.5, max_y * 1.1, '***', ha='center', fontsize=8)
```

## 使用不同的绘图库

### Matplotlib
- 对出版细节的最大控制权
- 最适合复杂的多面板图形
- 使用提供的样式文件来保持格式一致
- 有关详细示例，请参阅 `references/matplotlib_examples.md`

### 希博恩

Seaborn 为统计图形提供了一个高级的、面向数据集的接口，基于 matplotlib 构建。它擅长用最少的代码创建出版质量的统计可视化，同时保持与 matplotlib 自定义的完全兼容性。

**科学可视化的主要优势：**
- 自动统计估计和置信区间
- 内置支持多面板图形（分面）
- 默认情况下色盲友好的调色板
- 使用 pandas DataFrames 的面向数据集的 API
- 变量到视觉属性的语义映射

#### 出版风格快速入门

始终先应用 matplotlib 发布样式，然后配置 seaborn：

```python
import seaborn as sns
import matplotlib.pyplot as plt
from style_presets import apply_publication_style

# Apply publication style
apply_publication_style('default')

# Configure seaborn for publication
sns.set_theme(style='ticks', context='paper', font_scale=1.1)
sns.set_palette('colorblind')  # Use colorblind-safe palette

# Create figure
fig, ax = plt.subplots(figsize=(3.5, 2.5))
sns.scatterplot(data=df, x='time', y='response', 
                hue='treatment', style='condition', ax=ax)
sns.despine()  # Remove top and right spines
```

#### 出版物的常见绘图类型

**统计比较：**
```python
# Box plot with individual points for transparency
fig, ax = plt.subplots(figsize=(3.5, 3))
sns.boxplot(data=df, x='treatment', y='response', 
            order=['Control', 'Low', 'High'], palette='Set2', ax=ax)
sns.stripplot(data=df, x='treatment', y='response',
              order=['Control', 'Low', 'High'], 
              color='black', alpha=0.3, size=3, ax=ax)
ax.set_ylabel('Response (μM)')
sns.despine()
```

**分布分析：**
```python
# Violin plot with split comparison
fig, ax = plt.subplots(figsize=(4, 3))
sns.violinplot(data=df, x='timepoint', y='expression',
               hue='treatment', split=True, inner='quartile', ax=ax)
ax.set_ylabel('Gene Expression (AU)')
sns.despine()
```

**相关矩阵：**
```python
# Heatmap with proper colormap and annotations
fig, ax = plt.subplots(figsize=(5, 4))
corr = df.corr()
mask = np.triu(np.ones_like(corr, dtype=bool))  # Show only lower triangle
sns.heatmap(corr, mask=mask, annot=True, fmt='.2f',
            cmap='RdBu_r', center=0, square=True,
            linewidths=1, cbar_kws={'shrink': 0.8}, ax=ax)
plt.tight_layout()
```

**具有置信带的时间序列：**
```python
# Line plot with automatic CI calculation
fig, ax = plt.subplots(figsize=(5, 3))
sns.lineplot(data=timeseries, x='time', y='measurement',
             hue='treatment', style='replicate',
             errorbar=('ci', 95), markers=True, dashes=False, ax=ax)
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Measurement (AU)')
sns.despine()
```

#### Seaborn 的多面板人物

**使用 FacetGrid 进行自动分面：**
```python
# Create faceted plot
g = sns.relplot(data=df, x='dose', y='response',
                hue='treatment', col='cell_line', row='timepoint',
                kind='line', height=2.5, aspect=1.2,
                errorbar=('ci', 95), markers=True)
g.set_axis_labels('Dose (μM)', 'Response (AU)')
g.set_titles('{row_name} | {col_name}')
sns.despine()

# Save with correct DPI
from figure_export import save_publication_figure
save_publication_figure(g.figure, 'figure_facets', 
                       formats=['pdf', 'png'], dpi=300)
```

**将seaborn与matplotlib子图相结合：**
```python
# Create custom multi-panel layout
fig, axes = plt.subplots(2, 2, figsize=(7, 6))

# Panel A: Scatter with regression
sns.regplot(data=df, x='predictor', y='response', ax=axes[0, 0])
axes[0, 0].text(-0.15, 1.05, 'A', transform=axes[0, 0].transAxes,
                fontsize=10, fontweight='bold')

# Panel B: Distribution comparison
sns.violinplot(data=df, x='group', y='value', ax=axes[0, 1])
axes[0, 1].text(-0.15, 1.05, 'B', transform=axes[0, 1].transAxes,
                fontsize=10, fontweight='bold')

# Panel C: Heatmap
sns.heatmap(correlation_data, cmap='viridis', ax=axes[1, 0])
axes[1, 0].text(-0.15, 1.05, 'C', transform=axes[1, 0].transAxes,
                fontsize=10, fontweight='bold')

# Panel D: Time series
sns.lineplot(data=timeseries, x='time', y='signal', 
             hue='condition', ax=axes[1, 1])
axes[1, 1].text(-0.15, 1.05, 'D', transform=axes[1, 1].transAxes,
                fontsize=10, fontweight='bold')

plt.tight_layout()
sns.despine()
```

#### 出版物的调色板

Seaborn 包括几个色盲安全调色板：

```python
# Use built-in colorblind palette (recommended)
sns.set_palette('colorblind')

# Or specify custom colorblind-safe colors (Okabe-Ito)
okabe_ito = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000']
sns.set_palette(okabe_ito)

# For heatmaps and continuous data
sns.heatmap(data, cmap='viridis')  # Perceptually uniform
sns.heatmap(corr, cmap='RdBu_r', center=0)  # Diverging, centered
```

#### 在轴级函数和图形级函数之间进行选择

**轴级函数**（例如，`scatterplot`、`boxplot`、`heatmap`）：
- 在构建自定义多面板布局时使用
- 接受`ax=`参数以实现精确放置
- 与 matplotlib 子图更好地集成
- 更好地控制人物构图

```python
fig, ax = plt.subplots(figsize=(3.5, 2.5))
sns.scatterplot(data=df, x='x', y='y', hue='group', ax=ax)
```

**图形级函数**（例如，`relplot`、`catplot`、`displot`）：
- 用于按分类变量自动分面
- 创建具有一致风格的完整人物
- 非常适合探索性分析
- 使用 `height` 和 `aspect` 调整大小

```python
g = sns.relplot(data=df, x='x', y='y', col='category', kind='scatter')
```

#### Seaborn 的统计严谨性

Seaborn 自动计算并显示不确定性：

```python
# Line plot: shows mean ± 95% CI by default
sns.lineplot(data=df, x='time', y='value', hue='treatment',
             errorbar=('ci', 95))  # Can change to 'sd', 'se', etc.

# Bar plot: shows mean with bootstrapped CI
sns.barplot(data=df, x='treatment', y='response',
            errorbar=('ci', 95), capsize=0.1)

# Always specify error type in figure caption:
# "Error bars represent 95% confidence intervals"
```

#### 可供出版的 Seaborn 人物的最佳实践

1. **始终首先设置发布主题：**
   ```python
   sns.set_theme(style='ticks', context='paper', font_scale=1.1)
   ```

2. **使用色盲安全调色板：**
   ```python
   sns.set_palette('colorblind')
   ```

3. **删除不必要的元素：**
   ```python
   sns.despine()  # Remove top and right spines
   ```

4. **适当控制图形尺寸：**
   ```python
   # Axes-level: use matplotlib figsize
   fig, ax = plt.subplots(figsize=(3.5, 2.5))
   
   # Figure-level: use height and aspect
   g = sns.relplot(..., height=3, aspect=1.2)
   ```

5. **尽可能显示各个数据点：**
   ```python
   sns.boxplot(...)  # Summary statistics
   sns.stripplot(..., alpha=0.3)  # Individual points
   ```

6. **包含带有单位的正确标签：**
```python
   ax.set_xlabel('Time (hours)')
   ax.set_ylabel('Expression (AU)')
   ```

7. **以正确的分辨率导出：**
   ```python
   from figure_export import save_publication_figure
   save_publication_figure(fig, 'figure_name', 
                          formats=['pdf', 'png'], dpi=300)
   ```

#### 先进的 Seaborn 技术

**用于探索性分析的成对关系：**
```python
# Quick overview of all relationships
g = sns.pairplot(data=df, hue='condition', 
                 vars=['gene1', 'gene2', 'gene3'],
                 corner=True, diag_kind='kde', height=2)
```

**层次聚类热图：**
```python
# Cluster samples and features
g = sns.clustermap(expression_data, method='ward', 
                   metric='euclidean', z_score=0,
                   cmap='RdBu_r', center=0, 
                   figsize=(10, 8), 
                   row_colors=condition_colors,
                   cbar_kws={'label': 'Z-score'})
```

**带边际的联合分布：**
```python
# Bivariate distribution with context
g = sns.jointplot(data=df, x='gene1', y='gene2',
                  hue='treatment', kind='scatter',
                  height=6, ratio=4, marginal_kws={'kde': True})
```

#### Seaborn 常见问题和解决方案

**问题：情节区域外的图例**
```python
g = sns.relplot(...)
g._legend.set_bbox_to_anchor((0.9, 0.5))
```

**问题：标签重叠**
```python
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
```

**问题：最终尺寸的文字太小**
```python
sns.set_context('paper', font_scale=1.2)  # Increase if needed
```

#### 其他资源

有关 Seaborn 的更详细信息，请参阅：
- `scientific-packages/seaborn/SKILL.md` - 综合 Seaborn 文档
- `scientific-packages/seaborn/references/examples.md` - 实际用例
- `scientific-packages/seaborn/references/function_reference.md` - 完整的 API 参考
- `scientific-packages/seaborn/references/objects_interface.md` - 现代声明式 API

### 情节
- 用于探索的互动人物
- 导出静态图像以供发布
- 配置发布质量：
```python
fig.update_layout(
    font=dict(family='Arial, sans-serif', size=10),
    plot_bgcolor='white',
    # ... see matplotlib_examples.md Example 8
)
fig.write_image('figure.png', scale=3)  # scale=3 gives ~300 DPI
```

## 资源

### 参考目录

**根据需要加载这些以获取详细信息：**

- **`publication_guidelines.md`**：综合最佳实践
  - 分辨率和文件格式要求
  - 版式指南
  - 布局和构图规则
  - 统计严谨性要求
  - 完整的出版清单

- **`color_palettes.md`**：颜色使用指南
  - 具有 RGB 值的色盲友好调色板规范
  - 顺序和发散的色彩图建议
  - 可访问性测试程序
  - 特定领域的调色板（基因组学、显微镜学）

- **`journal_requirements.md`**：期刊特定规范
  - 出版商的技术要求
  - 文件格式和DPI规格
  - 人物尺寸要求
  - 快速参考表

- **`matplotlib_examples.md`**：实用代码示例
  - 10 个完整的工作示例
  - 线图、条形图、热图、多面板图
  - 期刊特定的图表示例
  - 每个库的提示（matplotlib、seaborn、plotly）

### 脚本目录

**使用这些帮助程序脚本进行自动化：**

- **`figure_export.py`**：导出实用程序
  - `save_publication_figure()`：以正确的 DPI 保存为多种格式
  - `save_for_journal()`：自动使用期刊特定要求
  - `check_figure_size()`：验证尺寸是否符合轴颈规格
  - 直接运行：`python scripts/figure_export.py` 为例

- **`style_presets.py`**：预配置样式
  - `apply_publication_style()`：应用预设样式（默认、自然、科学、单元格）
  - `set_color_palette()`：快速调色板切换
  - `configure_for_journal()`：单命令日志配置
  - 直接运行：`python scripts/style_presets.py`查看示例

### 资产目录

**在图中使用这些文件：**

- **`color_palettes.py`**：可导入的颜色定义
  - 所有推荐的调色板作为Python常量
  - `apply_palette()` 辅助函数
  - 可以直接导入到笔记本/脚本中

- **Matplotlib 样式文件**：与 `plt.style.use()` 一起使用
  - `publication.mplstyle`：一般发布质量
  - `nature.mplstyle`：自然期刊规范
  - `presentation.mplstyle`：海报/幻灯片的较大字体

## 工作流程总结

**创建出版物图的推荐工作流程：**

1. **计划**：确定目标期刊、图表类型和内容
2. **配置**：为日志应用适当的样式
   ```python
   from style_presets import configure_for_journal
   configure_for_journal('nature', 'single')
   ```
3. **创建**：使用适当的标签、颜色、统计数据构建图形
4. **验证**：检查尺寸、字体、颜色、可访问性
   ```python
   from figure_export import check_figure_size
   check_figure_size(fig, journal='nature')
   ```
5. **导出**：以所需格式保存
   ```python
   from figure_export import save_for_journal
   save_for_journal(fig, 'figure1', 'nature', 'combination')
   ```
6. **审查**：在手稿上下文中查看最终尺寸

## 要避免的常见陷阱

1. **字体太小**：以最终尺寸打印时文本无法阅读
2. **JPEG 格式**：切勿将 JPEG 用于图形/绘图（创建伪影）
3. **红绿颜色**：约8%的男性无法区分
4. **低分辨率**：出版物中的像素化数字
5. **缺少单位**：始终用单位标记轴
6. **3D效果**：扭曲感知，完全避免
7. **图表垃圾**：删除不必要的网格线、装饰
8. **截断轴**：条形图从零开始，除非有科学依据
9. **样式不一致**：同一手稿中的各个图形的字体/颜色不同
10. **无误差线**：始终表现出不确定性

## 最终清单

在提交数据之前，请验证：

- [ ] 分辨率符合期刊要求（300+ DPI）
- [ ] 文件格式正确（绘图为矢量，图像为 TIFF）
- [ ] 图形尺寸符合期刊规格
- [ ] 所有文本均以最终尺寸可读（≥6 pt）
- [ ] 颜色对色盲者友好
- [ ] 图形以灰度显示
- [ ] 所有标有单位的轴
- [ ] 错误栏与标题中的定义一起出现
- [ ] 面板标签存在且一致
- [ ] 没有图表垃圾或 3D 效果
- [ ] 所有图形的字体一致
- [ ] 明确标记统计显着性
- [ ]图例清晰完整

使用此技能可确保科学数据符合最高的出版标准，同时保持所有读者都可以访问。