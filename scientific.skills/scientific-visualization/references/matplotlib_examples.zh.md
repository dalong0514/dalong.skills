<!-- 此文件由机器翻译自 matplotlib_examples.md -->

# 可供发布的 Matplotlib 示例

## 概述

本参考提供了使用 Matplotlib、Seaborn 和 Plotly 创建可发表的科学图形的实用代码示例。所有示例均遵循 `publication_guidelines.md` 中的最佳实践，并使用 `color_palettes.md` 中的色盲友好调色板。

## 设置和配置

### 出版质量的 Matplotlib 配置

```python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# Set publication quality parameters
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['font.size'] = 8
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']
mpl.rcParams['axes.labelsize'] = 9
mpl.rcParams['axes.titlesize'] = 9
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
mpl.rcParams['legend.fontsize'] = 7
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['lines.linewidth'] = 1.5

# Use colorblind-friendly colors (Okabe-Ito palette)
okabe_ito = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000']
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=okabe_ito)

# Use perceptually uniform colormap
mpl.rcParams['image.cmap'] = 'viridis'
```

### 保存辅助函数

<<<代码块_1>>>

## 示例 1：带误差线的线图

<<<代码块_2>>>

## 示例 2：多面板图

<<<代码块_3>>>

## 示例 3：包含各个点的箱线图

<<<代码块_4>>>

## 示例 4：带颜色条的热图

<<<代码块_5>>>

## 示例 5：Seaborn 小提琴图

<<<代码块_6>>>

## 示例 6：带有回归的科学散点图

```python
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# Generate data with correlation
np.random.seed(42)
x = np.random.randn(100)
y = 2.5 * x + np.random.randn(100) * 0.8

# Calculate regression
slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

# Create figure
fig, ax = plt.subplots(figsize=(3.5, 3.5))

# Scatter plot
ax.scatter(x, y, s=15, alpha=0.6, color='#0072B2', edgecolors='none')

# Regression line
x_line = np.array([x.min(), x.max()])
y_line = slope * x_line + intercept
ax.plot(x_line, y_line, 'r-', linewidth=1.5, label=f'y = {slope:.2f}x + {intercept:.2f}')

# Add statistics text
stats_text = f'$R^2$ = {r_value**2:.3f}\n$p$ < 0.001' if p_value < 0.001 else f'$R^2$ = {r_value**2:.3f}\n$p$ = {p_value:.3f}'
ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
        verticalalignment='top', fontsize=7,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5))

# Customize
ax.set_xlabel('Predictor variable')
ax.set_ylabel('Response variable')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

fig.tight_layout()
save_publication_figure(fig, 'scatter_regression')
plt.show()
```

## 示例 7：具有阴影误差的时间序列

```python
import matplotlib.pyplot as plt
import numpy as np

# Generate time series data
np.random.seed(42)
time = np.linspace(0, 24, 100)
n_replicates = 5

# Simulate multiple replicates
data = np.array([10 * np.exp(-time/10) + np.random.normal(0, 0.5, 100)
                 for _ in range(n_replicates)])

# Calculate mean and SEM
mean = data.mean(axis=0)
sem = data.std(axis=0) / np.sqrt(n_replicates)

# Create figure
fig, ax = plt.subplots(figsize=(4, 2.5))

# Plot mean line
ax.plot(time, mean, linewidth=1.5, color='#0072B2', label='Mean ± SEM')

# Add shaded error region
ax.fill_between(time, mean - sem, mean + sem,
                alpha=0.3, color='#0072B2', linewidth=0)

# Customize
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Concentration (μM)')
ax.legend(frameon=False, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(0, 24)
ax.set_ylim(0, 12)

fig.tight_layout()
save_publication_figure(fig, 'timeseries_shaded')
plt.show()
```

## 示例 8：Plotly 交互式图形

```python
import plotly.graph_objects as go
import numpy as np

# Generate data
np.random.seed(42)
x = np.random.randn(100)
y = 2*x + np.random.randn(100)
colors = np.random.choice(['Group A', 'Group B'], 100)

# Okabe-Ito colors for Plotly
okabe_ito_plotly = ['#E69F00', '#56B4E9']

# Create figure
fig = go.Figure()

for group, color in zip(['Group A', 'Group B'], okabe_ito_plotly):
    mask = colors == group
    fig.add_trace(go.Scatter(
        x=x[mask], y=y[mask],
        mode='markers',
        name=group,
        marker=dict(size=6, color=color, opacity=0.6)
    ))

# Update layout for publication quality
fig.update_layout(
    width=500,
    height=400,
    font=dict(family='Arial, sans-serif', size=10),
    plot_bgcolor='white',
    xaxis=dict(
        title='Variable X',
        showgrid=False,
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=False
    ),
    yaxis=dict(
        title='Variable Y',
        showgrid=False,
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=False
    ),
    legend=dict(
        x=0.02,
        y=0.98,
        bgcolor='rgba(255,255,255,0.8)',
        bordercolor='gray',
        borderwidth=0.5
    )
)

# Save as static image (requires kaleido)
fig.write_image('plotly_scatter.png', width=500, height=400, scale=3)  # scale=3 gives ~300 DPI
fig.write_html('plotly_scatter.html')  # Interactive version

fig.show()
```

## 示例 9：具有显着性的分组条形图

```python
import matplotlib.pyplot as plt
import numpy as np

# Data
categories = ['WT', 'Mutant A', 'Mutant B']
control_means = [100, 85, 70]
control_sem = [5, 6, 5]
treatment_means = [100, 120, 140]
treatment_sem = [6, 8, 9]

x = np.arange(len(categories))
width = 0.35

fig, ax = plt.subplots(figsize=(3.5, 3))

# Create bars
bars1 = ax.bar(x - width/2, control_means, width, yerr=control_sem,
               capsize=3, label='Control', color='#0072B2', alpha=0.8)
bars2 = ax.bar(x + width/2, treatment_means, width, yerr=treatment_sem,
               capsize=3, label='Treatment', color='#E69F00', alpha=0.8)

# Add significance markers
def add_significance_bar(ax, x1, x2, y, h, text):
    """Add significance bar between two bars"""
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], linewidth=0.8, c='black')
    ax.text((x1+x2)/2, y+h, text, ha='center', va='bottom', fontsize=7)

# Mark significant differences
add_significance_bar(ax, x[1]-width/2, x[1]+width/2, 135, 3, '***')
add_significance_bar(ax, x[2]-width/2, x[2]+width/2, 155, 3, '***')

# Customize
ax.set_ylabel('Activity (% of WT control)')
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend(frameon=False, loc='upper left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0, 180)

# Add note about significance
ax.text(0.98, 0.02, '*** p < 0.001', transform=ax.transAxes,
        ha='right', va='bottom', fontsize=6)

fig.tight_layout()
save_publication_figure(fig, 'grouped_bar_significance')
plt.show()
```

## 示例 10：可供发表的 Nature 图

```python
import matplotlib.pyplot as plt
import numpy as np
from string import ascii_lowercase

# Nature specifications: 89mm single column
inch_per_mm = 0.0393701
width_mm = 89
height_mm = 110
figsize = (width_mm * inch_per_mm, height_mm * inch_per_mm)

fig = plt.figure(figsize=figsize)
gs = fig.add_gridspec(3, 2, hspace=0.5, wspace=0.4,
                      left=0.12, right=0.95, top=0.96, bottom=0.08)

# Panel a: Time course
ax_a = fig.add_subplot(gs[0, :])
time = np.linspace(0, 48, 100)
for i, label in enumerate(['Control', 'Treatment']):
    y = (1 + i*0.5) * np.exp(-time/20) * (1 + 0.3*np.sin(time/5))
    ax_a.plot(time, y, linewidth=1.2, label=label)
ax_a.set_xlabel('Time (h)', fontsize=7)
ax_a.set_ylabel('Growth (OD$_{600}$)', fontsize=7)
ax_a.legend(frameon=False, fontsize=6)
ax_a.tick_params(labelsize=6)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# Panel b: Bar plot
ax_b = fig.add_subplot(gs[1, 0])
categories = ['A', 'B', 'C']
values = [1.0, 1.5, 2.2]
errors = [0.1, 0.15, 0.2]
ax_b.bar(categories, values, yerr=errors, capsize=2, width=0.6,
         color='#0072B2', alpha=0.8)
ax_b.set_ylabel('Fold change', fontsize=7)
ax_b.tick_params(labelsize=6)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# Panel c: Heatmap
ax_c = fig.add_subplot(gs[1, 1])
data = np.random.randn(8, 6)
im = ax_c.imshow(data, cmap='viridis', aspect='auto')
ax_c.set_xlabel('Sample', fontsize=7)
ax_c.set_ylabel('Gene', fontsize=7)
ax_c.tick_params(labelsize=6)

# Panel d: Scatter
ax_d = fig.add_subplot(gs[2, :])
x = np.random.randn(50)
y = 2*x + np.random.randn(50)*0.5
ax_d.scatter(x, y, s=8, alpha=0.6, color='#E69F00')
ax_d.set_xlabel('Expression gene X', fontsize=7)
ax_d.set_ylabel('Expression gene Y', fontsize=7)
ax_d.tick_params(labelsize=6)
ax_d.spines['top'].set_visible(False)
ax_d.spines['right'].set_visible(False)

# Add lowercase panel labels (Nature style)
for i, ax in enumerate([ax_a, ax_b, ax_c, ax_d]):
    ax.text(-0.2, 1.1, f'{ascii_lowercase[i]}', transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top')

# Save in Nature-preferred format
fig.savefig('nature_figure.pdf', dpi=1000, bbox_inches='tight',
           facecolor='white', edgecolor='none')
fig.savefig('nature_figure.png', dpi=300, bbox_inches='tight',
           facecolor='white', edgecolor='none')

plt.show()
```

## 每个图书馆的提示

### Matplotlib
- 使用`fig.tight_layout()`或`constrained_layout=True`来防止重叠
- 将 DPI 设置为 300-600 进行发布
- 使用矢量格式（PDF、EPS）绘制线图
- 在 PDF/EPS 文件中嵌入字体

### 希博恩
- 基于 matplotlib 构建，因此所有 matplotlib 自定义都可以工作
- 使用 `sns.set_style('ticks')` 或 `'whitegrid'` 获得干净的外观
- `sns.despine()` 删除顶部和右侧的书脊
- 使用 `sns.set_palette()` 设置自定义调色板

### 情节
- 非常适合交互式探索性分析
- 使用`fig.write_image()`导出静态图像（需要kaleido包）
- 使用`scale`参数来控制DPI（scale=3 ≈ 300 DPI）
- 广泛更新布局以提高出版物质量

## 通用工作流程

1. **使用默认设置探索**
2. **应用发布配置**（请参阅设置部分）
3. **创建适当大小的图**（检查期刊要求）
4. **自定义颜色**（使用色盲友好的调色板）
5. **调整字体和线宽**（以最终尺寸可读）
6. **删除图表垃圾**（顶部/右侧的脊柱，过多的网格）
7. **添加带有单位的清晰标签**
8. **灰度测试**
9. **以多种格式保存**（矢量为 PDF，光栅为 PNG）
10. **在最终上下文中验证**（导入手稿以检查尺寸）

## 资源

- Matplotlib 文档：https://matplotlib.org/
- Seaborn 画廊：https://seaborn.pydata.org/examples/index.html
- Plotly 文档：https://plotly.com/python/
- NatureMethods Points of View：数据可视化专栏存档