<!-- 此文件由机器翻译自 SKILL.md -->

---
姓名：海博恩
描述：“统计可视化。散点图、箱线图、小提琴图、热图、配对图、回归、相关矩阵、KDE、多面图，用于探索性分析和出版数据。”
---

# Seaborn 统计可视化

## 概述

Seaborn 是一个 Python 可视化库，用于创建出版质量的统计图形。使用此技能以最少的代码进行面向数据集的绘图、多变量分析、自动统计估计和复杂的多面板图形。

## 设计理念

Seaborn 遵循以下核心原则：

1. **面向数据集**：直接使用DataFrames和命名变量而不是抽象坐标
2. **语义映射**：自动将数据值转换为视觉属性（颜色、大小、样式）
3. **统计意识**：内置聚合、误差估计和置信区间
4. **审美默认设置**：开箱即用的出版就绪主题和调色板
5. **Matplotlib集成**：在需要时与matplotlib定制完全兼容

## 快速入门

```python
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Load example dataset
df = sns.load_dataset('tips')

# Create a simple visualization
sns.scatterplot(data=df, x='total_bill', y='tip', hue='day')
plt.show()
```

## 核心绘图接口

### 函数接口（传统）

函数接口提供了按可视化类型组织的专门绘图函数。每个类别都有**轴级**函数（绘制到单轴）和**图形级**函数（通过分面管理整个图形）。

**何时使用：**
- 快速探索性分析
- 单一目的可视化
- 当您需要特定的绘图类型时

### 对象接口（现代）

`seaborn.objects` 接口提供了类似于 ggplot2 的声明性、可组合 API。通过链接方法来指定数据映射、标记、转换和比例来构建可视化。

**何时使用：**
- 复杂的分层可视化
- 当您需要对转换进行细粒度控制时
- 构建自定义绘图类型
- 程序化绘图生成

<<<代码块_1>>>

## 按类别绘制函数

### 关系图（变量之间的关系）

**用途：**探索两个或多个变量如何相互关联

- `scatterplot()` - 将单个观察结果显示为点
- `lineplot()` - 显示趋势和变化（自动聚合和计算 CI）
- `relplot()` - 具有自动分面功能的图形级界面

**关键参数：**
- `x`, `y` - 主变量
- `hue` - 附加分类/连续变量的颜色编码
- `size` - 点/线尺寸编码
- `style` - 标记/线条样式编码
- `col`, `row` - 分为多个子图（仅限图形级别）

<<<代码块_2>>>

### 分布图（单变量分布和双变量分布）

**用于：** 了解数据分布、形状和概率密度

- `histplot()` - 基于条形的频率分布，具有灵活的分级功能
- `kdeplot()` - 使用高斯核平滑密度估计
- `ecdfplot()` - 经验累积分布（无需调整参数）
- `rugplot()` - 单独观察刻度线
- `displot()` - 单变量和双变量分布的图形级界面
- `jointplot()` - 具有边际分布的双变量图
- `pairplot()` - 跨数据集的成对关系矩阵

**关键参数：**
- `x`, `y` - 变量（对于单变量，y 可选）
- `hue` - 按类别单独分布
- `stat` - 标准化：“计数”、“频率”、“概率”、“密度”
- `bins` / `binwidth` - 直方图分级控制
- `bw_adjust` - KDE 带宽乘数（越高 = 越平滑）
- `fill` - 填充曲线下区域
- `multiple` - 如何处理色调：“layer”、“stack”、“dodge”、“fill”

<<<代码块_3>>>

### 分类图（跨类别比较）

**用于：** 比较离散类别的分布或统计数据

**分类散点图：**
- `stripplot()` - 带抖动的点显示所有观察结果
- `swarmplot()` - 非重叠点（beeswarm 算法）

**分布比较：**
- `boxplot()` - 四分位数和离群值
- `violinplot()` - KDE + 四分位数信息
- `boxenplot()` - 针对较大数据集的增强箱线图

**统计估计：**
- `barplot()` - 具有置信区间的平均值/合计
- `pointplot()` - 带连接线的点估计
- `countplot()` - 每个类别的观察计数

**人物级别：**
- `catplot()` - 分面分类图（设置 `kind` 参数）

**关键参数：**
- `x`, `y` - 变量（通常是分类变量）
- `hue` - 附加分类分组
- `order`, `hue_order` - 控制类别排序
- `dodge` - 并排单独的色调级别
- `orient` - “v”（垂直）或“h”（水平）
- `kind` - catplot 的绘图类型：“strip”、“swarm”、“box”、“violin”、“bar”、“point”

<<<代码块_4>>>

### 回归图（线性关系）

**用途：** 可视化线性回归和残差

- `regplot()` - 带有散点+拟合线的轴级回归图
- `lmplot()` - 具有分面支持的图形级别
- `residplot()` - 用于评估模型拟合的残差图

**关键参数：**
- `x`, `y` - 要回归的变量
- `order` - 多项式回归阶数
- `logistic` - 拟合逻辑回归
- `robust` - 使用稳健回归（对异常值不太敏感）
- `ci` - 置信区间宽度（默认 95）
- `scatter_kws`、`line_kws` - 自定义散点和线条属性

<<<代码块_5>>>

### 矩阵图（矩形数据）

**用途：** 可视化矩阵、相关性和网格结构数据

- `heatmap()` - 带注释的颜色编码矩阵
- `clustermap()` - 分层集群热图

**关键参数：**
- `data` - 2D 矩形数据集（DataFrame 或数组）
- `annot` - 显示单元格中的值
- `fmt` - 注释的格式字符串（例如“.2f”）
- `cmap` - 颜色图名称
- `center` - 颜色图中心的值（用于发散颜色图）
- `vmin`、`vmax` - 色阶限制
- `square` - 强制方形单元格
- `linewidths` - 单元格之间的间隙

<<<代码块_6>>>

## 多图网格

Seaborn 提供了用于创建复杂的多面板图形的网格对象：

### FacetGrid

根据分类变量创建子图。通过图形级函数（`relplot`、`displot`、`catplot`）调用时最有用，但可直接用于自定义绘图。

```python
g = sns.FacetGrid(df, col='time', row='sex', hue='smoker')
g.map(sns.scatterplot, 'total_bill', 'tip')
g.add_legend()
```

### 配对网格

显示数据集中所有变量之间的成对关系。

```python
g = sns.PairGrid(df, hue='species')
g.map_upper(sns.scatterplot)
g.map_lower(sns.kdeplot)
g.map_diag(sns.histplot)
g.add_legend()
```

### 联合网格

将二变量图与边际分布相结合。

```python
g = sns.JointGrid(data=df, x='total_bill', y='tip')
g.plot_joint(sns.scatterplot)
g.plot_marginals(sns.histplot)
```

## 图形级函数与轴级函数

理解这种区别对于有效使用seaborn至关重要：

### 轴级函数
- 绘制到单个 matplotlib `Axes` 对象
- 轻松集成到复杂的 matplotlib 图形中
- 接受`ax=`参数以实现精确放置
- 返回`Axes`对象
- 示例：`scatterplot`、`histplot`、`boxplot`、`regplot`、`heatmap`

**何时使用：**
- 构建自定义多图布局
- 结合不同的绘图类型
- 需要matplotlib级别的控制
- 与现有的 matplotlib 代码集成

```python
fig, axes = plt.subplots(2, 2, figsize=(10, 10))
sns.scatterplot(data=df, x='x', y='y', ax=axes[0, 0])
sns.histplot(data=df, x='x', ax=axes[0, 1])
sns.boxplot(data=df, x='cat', y='y', ax=axes[1, 0])
sns.kdeplot(data=df, x='x', y='y', ax=axes[1, 1])
```

### 图形级函数
- 管理整个图形，包括所有子图
- 通过 `col` 和 `row` 参数进行内置分面
- 返回 `FacetGrid`、`JointGrid` 或 `PairGrid` 对象
- 使用 `height` 和 `aspect` 调整大小（每个子图）
- 无法放置在现有图形中
- 示例：`relplot`、`displot`、`catplot`、`lmplot`、`jointplot`、`pairplot`

**何时使用：**
- 多面可视化（小倍数）
- 快速探索性分析
- 一致的多面板布局
- 不需要与其他绘图类型结合

```python
# Automatic faceting
sns.relplot(data=df, x='x', y='y', col='category', row='group',
            hue='type', height=3, aspect=1.2)
```

## 数据结构要求

### 长格式数据（首选）

每个变量是一列，每个观察值是一行。这种“整洁”的格式提供了最大的灵活性：

```python
# Long-form structure
   subject  condition  measurement
0        1    control         10.5
1        1  treatment         12.3
2        2    control          9.8
3        2  treatment         13.1
```

**优点：**
- 适用于所有seaborn功能
- 轻松将变量重新映射到视觉属性
- 支持任意复杂度
- 自然的 DataFrame 操作

### 宽格式数据
变量分布在各列中。对于简单的矩形数据很有用：

```python
# Wide-form structure
   control  treatment
0     10.5       12.3
1      9.8       13.1
```

**使用案例：**
- 简单的时间序列
- 相关矩阵
- 热图
- 阵列数据的快速绘图

**将宽转换为长：**
```python
df_long = df.melt(var_name='condition', value_name='measurement')
```

## 调色板

Seaborn 为不同的数据类型提供精心设计的调色板：

### 定性调色板（分类数据）

通过色调变化区分类别：
- `"deep"` - 默认、鲜艳的颜色
- `"muted"` - 更柔和，饱和度更低
- `"pastel"` - 浅色，去饱和
- `"bright"` - 高度饱和
- `"dark"` - 暗值
- `"colorblind"` - 对于色觉缺陷来说是安全的

```python
sns.set_palette("colorblind")
sns.color_palette("Set2")
```

### 顺序调色板（有序数据）

显示从低值到高值的进展：
- `"rocket"`、`"mako"` - 宽亮度范围（适合热图）
- `"flare"`、`"crest"` - 亮度受限（适用于点/线）
- `"viridis"`、`"magma"`、`"plasma"` - Matplotlib 感知均匀

```python
sns.heatmap(data, cmap='rocket')
sns.kdeplot(data=df, x='x', y='y', cmap='mako', fill=True)
```

### 发散调色板（中心数据）

强调与中点的偏差：
- `"vlag"` - 蓝色到红色
- `"icefire"` - 蓝色到橙色
- `"coolwarm"` - 冷到暖
- `"Spectral"` - 彩虹发散

```python
sns.heatmap(correlation_matrix, cmap='vlag', center=0)
```

### 自定义调色板

```python
# Create custom palette
custom = sns.color_palette("husl", 8)

# Light to dark gradient
palette = sns.light_palette("seagreen", as_cmap=True)

# Diverging palette from hues
palette = sns.diverging_palette(250, 10, as_cmap=True)
```

## 主题和美学

### 设置主题

`set_theme()` 控制整体外观：

```python
# Set complete theme
sns.set_theme(style='whitegrid', palette='pastel', font='sans-serif')

# Reset to defaults
sns.set_theme()
```

### 风格

控制背景和网格外观：
- `"darkgrid"` - 带白色网格的灰色背景（默认）
- `"whitegrid"` - 带灰色网格的白色背景
- `"dark"` - 灰色背景，无网格
- `"white"` - 白色背景，无网格
- `"ticks"` - 带轴刻度的白色背景

```python
sns.set_style("whitegrid")

# Remove spines
sns.despine(left=False, bottom=False, offset=10, trim=True)

# Temporary style
with sns.axes_style("white"):
    sns.scatterplot(data=df, x='x', y='y')
```

### 上下文

针对不同用例缩放元素：
- `"paper"` - 最小（默认）
- `"notebook"` - 稍大
- `"talk"` - 演示幻灯片
- `"poster"` - 大格式

```python
sns.set_context("talk", font_scale=1.2)

# Temporary context
with sns.plotting_context("poster"):
    sns.barplot(data=df, x='category', y='value')
```

## 最佳实践

### 1. 数据准备

始终使用结构良好且具有有意义的列名称的 DataFrame：

```python
# Good: Named columns in DataFrame
df = pd.DataFrame({'bill': bills, 'tip': tips, 'day': days})
sns.scatterplot(data=df, x='bill', y='tip', hue='day')

# Avoid: Unnamed arrays
sns.scatterplot(x=x_array, y=y_array)  # Loses axis labels
```

### 2. 选择正确的绘图类型

**连续x、连续y：** `scatterplot`、`lineplot`、`kdeplot`、`regplot`
**连续 x，分类 y：** `violinplot`、`boxplot`、`stripplot`、`swarmplot`
**一个连续变量：** `histplot`、`kdeplot`、`ecdfplot`
**相关性/矩阵：** `heatmap`, `clustermap`
**成对关系：** `pairplot`、`jointplot`

### 3. 使用图形级函数进行分面

```python
# Instead of manual subplot creation
sns.relplot(data=df, x='x', y='y', col='category', col_wrap=3)

# Not: Creating subplots manually for simple faceting
```

### 4. 利用语义映射

使用 `hue`、`size` 和 `style` 对附加维度进行编码：

```python
sns.scatterplot(data=df, x='x', y='y',
                hue='category',      # Color by category
                size='importance',    # Size by continuous variable
                style='type')         # Marker style by type
```

### 5. 控制统计估计

许多函数自动计算统计数据。了解并定制：

```python
# Lineplot computes mean and 95% CI by default
sns.lineplot(data=df, x='time', y='value',
             errorbar='sd')  # Use standard deviation instead

# Barplot computes mean by default
sns.barplot(data=df, x='category', y='value',
            estimator='median',  # Use median instead
            errorbar=('ci', 95))  # Bootstrapped CI
```

### 6.与Matplotlib结合

Seaborn 与 matplotlib 无缝集成以进行微调：

```python
ax = sns.scatterplot(data=df, x='x', y='y')
ax.set(xlabel='Custom X Label', ylabel='Custom Y Label',
       title='Custom Title')
ax.axhline(y=0, color='r', linestyle='--')
plt.tight_layout()
```

### 7. 保存高质量的数据

```python
fig = sns.relplot(data=df, x='x', y='y', col='group')
fig.savefig('figure.png', dpi=300, bbox_inches='tight')
fig.savefig('figure.pdf')  # Vector format for publications
```

## 常见模式

### 探索性数据分析

```python
# Quick overview of all relationships
sns.pairplot(data=df, hue='target', corner=True)

# Distribution exploration
sns.displot(data=df, x='variable', hue='group',
            kind='kde', fill=True, col='category')

# Correlation analysis
corr = df.corr()
sns.heatmap(corr, annot=True, cmap='coolwarm', center=0)
```

### 出版质量数据

```python
sns.set_theme(style='ticks', context='paper', font_scale=1.1)

g = sns.catplot(data=df, x='treatment', y='response',
                col='cell_line', kind='box', height=3, aspect=1.2)
g.set_axis_labels('Treatment Condition', 'Response (μM)')
g.set_titles('{col_name}')
sns.despine(trim=True)

g.savefig('figure.pdf', dpi=300, bbox_inches='tight')
```

### 复杂的多面板人物

```python
# Using matplotlib subplots with seaborn
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

sns.scatterplot(data=df, x='x1', y='y', hue='group', ax=axes[0, 0])
sns.histplot(data=df, x='x1', hue='group', ax=axes[0, 1])
sns.violinplot(data=df, x='group', y='y', ax=axes[1, 0])
sns.heatmap(df.pivot_table(values='y', index='x1', columns='x2'),
            ax=axes[1, 1], cmap='viridis')

plt.tight_layout()
```

### 具有置信带的时间序列

```python
# Lineplot automatically aggregates and shows CI
sns.lineplot(data=timeseries, x='date', y='measurement',
             hue='sensor', style='location', errorbar='sd')

# For more control
g = sns.relplot(data=timeseries, x='date', y='measurement',
                col='location', hue='sensor', kind='line',
                height=4, aspect=1.5, errorbar=('ci', 95))
g.set_axis_labels('Date', 'Measurement (units)')
```

## 故障排除

### 问题：绘图区域外的图例

默认情况下，图形级函数将图例放置在外部。移入内部：

```python
g = sns.relplot(data=df, x='x', y='y', hue='category')
g._legend.set_bbox_to_anchor((0.9, 0.5))  # Adjust position
```

### 问题：重叠标签

```python
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
```

### 问题：图形太小

对于图形级函数：
```python
sns.relplot(data=df, x='x', y='y', height=6, aspect=1.5)
```

对于轴级函数：
```python
fig, ax = plt.subplots(figsize=(10, 6))
sns.scatterplot(data=df, x='x', y='y', ax=ax)
```

### 问题：颜色不够鲜明

```python
# Use a different palette
sns.set_palette("bright")

# Or specify number of colors
palette = sns.color_palette("husl", n_colors=len(df['category'].unique()))
sns.scatterplot(data=df, x='x', y='y', hue='category', palette=palette)
```

### 问题：KDE 太平滑或锯齿状

```python
# Adjust bandwidth
sns.kdeplot(data=df, x='x', bw_adjust=0.5)  # Less smooth
sns.kdeplot(data=df, x='x', bw_adjust=2)    # More smooth
```

## 资源

该技能包括用于更深入探索的参考材料：

###参考资料/

- `function_reference.md` - 所有seaborn函数的综合列表以及参数和示例
- `objects_interface.md` - 现代 seaborn.objects API 的详细指南
- `examples.md` - 不同分析场景的常见用例和代码模式
根据需要加载详细函数签名、高级参数或特定示例的参考文件。