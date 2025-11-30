<!-- 此文件由机器翻译自 styling_guide.md -->

# Matplotlib 样式指南

有关样式化和自定义 matplotlib 可视化的综合指南。

## 色彩图

### 色彩图类别

**1.感知均匀序列**
最适合从低值到高值的有序数据。
- `viridis`（默认，色盲友好）
- `plasma`
- `inferno`
- `magma`
- `cividis`（针对色盲观众进行了优化）

**用途：**
```python
im = ax.imshow(data, cmap='viridis')
scatter = ax.scatter(x, y, c=values, cmap='plasma')
```

**2.顺序**
用于有序数据的传统颜色图。
- `Blues`、`Greens`、`Reds`、`Oranges`、`Purples`
- `YlOrBr`、`YlOrRd`、`OrRd`、`PuRd`
- `BuPu`、`GnBu`、`PuBu`、`YlGnBu`

**3.发散**
最适合具有有意义的中心点（例如零、平均值）的数据。
- `coolwarm`（蓝色到红色）
- `RdBu`（红-蓝）
- `RdYlBu`（红-黄-蓝）
- `RdYlGn`（红-黄-绿）
- `PiYG`、`PRGn`、`BrBG`、`PuOr`、`RdGy`

**用途：**
<<<代码块_1>>>

**4.定性**
最适合没有固有排序的分类/名义数据。
- `tab10`（10 种不同的颜色）
- `tab20`（20 种不同的颜色）
- `Set1`、`Set2`、`Set3`
- `Pastel1`、`Pastel2`
- `Dark2`、`Accent`、`Paired`

**用途：**
<<<代码块_2>>>

**5.循环**
最适合循环数据（例如相位、角度）。
- `twilight`
- `twilight_shifted`
- `hsv`

### 色彩图最佳实践

1. **避免 `jet` 颜色图** - 感知上不统一，具有误导性
2. **使用感知统一的颜色图** - `viridis`、`plasma`、`cividis`
3. **考虑色盲用户** - 使用 `viridis`、`cividis`，或使用色盲模拟器进行测试
4. **将颜色图与数据类型匹配**：
   - 顺序：增加/减少数据
   - 发散：具有有意义中心的数据
   - 定性：类别
5. **反向颜色图** - 添加 `_r` 后缀：`viridis_r`、`coolwarm_r`

### 创建自定义颜色图

<<<代码块_3>>>

### 离散色彩图

<<<代码块_4>>>

## 样式表

### 使用内置样式

<<<代码块_5>>>

### 流行的内置样式

- `default` - Matplotlib 的默认样式
- `classic` - 经典 matplotlib 外观（2.0 之前）
- `seaborn-v0_8-*` - Seaborn 风格
  - `seaborn-v0_8-darkgrid`、`seaborn-v0_8-whitegrid`
  - `seaborn-v0_8-dark`、`seaborn-v0_8-white`
  - `seaborn-v0_8-ticks`、`seaborn-v0_8-poster`、`seaborn-v0_8-talk`
- `ggplot` - ggplot2 风格
- `bmh` - 黑客风格的贝叶斯方法
- `fivethirtyeight` - FiveThirtyEight 样式
- `grayscale` - 灰度样式

### 创建自定义样式表

创建一个名为 `custom_style.mplstyle` 的文件：

<<<代码块_6>>>

加载及使用：
```python
plt.style.use('path/to/custom_style.mplstyle')
```

## rcParams 配置

### 全局配置

```python
import matplotlib.pyplot as plt

# Configure globally
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14

# Or update multiple at once
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'lines.linewidth': 2
})
```

### 临时配置

```python
# Context manager for temporary changes
with plt.rc_context({'font.size': 14, 'lines.linewidth': 2.5}):
    fig, ax = plt.subplots()
    ax.plot(x, y)
```

### 常用 rcParams

**图形设置：**
```python
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['figure.dpi'] = 100
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.edgecolor'] = 'white'
plt.rcParams['figure.autolayout'] = False
plt.rcParams['figure.constrained_layout.use'] = True
```

**字体设置：**
```python
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.size'] = 12
plt.rcParams['font.weight'] = 'normal'
```

**轴设置：**
```python
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelweight'] = 'normal'
plt.rcParams['axes.spines.top'] = True
plt.rcParams['axes.spines.right'] = True
```

**线路设置：**
```python
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.linestyle'] = '-'
plt.rcParams['lines.marker'] = 'None'
plt.rcParams['lines.markersize'] = 6
```

**保存设置：**
```python
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.format'] = 'png'
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['savefig.pad_inches'] = 0.1
plt.rcParams['savefig.transparent'] = False
```

## 调色板

### 命名颜色集

```python
# Tableau colors
tableau_colors = plt.cm.tab10.colors

# CSS4 colors (subset)
css_colors = ['steelblue', 'coral', 'teal', 'goldenrod', 'crimson']

# Manual definition
custom_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
```

### 颜色循环

```python
# Set default color cycle
from cycler import cycler
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
plt.rcParams['axes.prop_cycle'] = cycler(color=colors)

# Or combine color and line style
plt.rcParams['axes.prop_cycle'] = cycler(color=colors) + cycler(linestyle=['-', '--', ':', '-.'])
```

### 调色板生成

```python
# Evenly spaced colors from colormap
n_colors = 5
colors = plt.cm.viridis(np.linspace(0, 1, n_colors))

# Use in plot
for i, (x, y) in enumerate(data):
    ax.plot(x, y, color=colors[i])
```

## 版式

### 字体配置

```python
# Set font family
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']

# Or sans-serif
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']

# Or monospace
plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.monospace'] = ['Courier New', 'DejaVu Sans Mono']
```

### 文本中的字体属性

```python
from matplotlib import font_manager

# Specify font properties
ax.text(x, y, 'Text',
        fontsize=14,
        fontweight='bold',  # 'normal', 'bold', 'heavy', 'light'
        fontstyle='italic',  # 'normal', 'italic', 'oblique'
        fontfamily='serif')

# Use specific font file
prop = font_manager.FontProperties(fname='path/to/font.ttf')
ax.text(x, y, 'Text', fontproperties=prop)
```

### 数学文本

```python
# LaTeX-style math
ax.set_title(r'$\alpha > \beta$')
ax.set_xlabel(r'$\mu \pm \sigma$')
ax.text(x, y, r'$\int_0^\infty e^{-x} dx = 1$')

# Subscripts and superscripts
ax.set_ylabel(r'$y = x^2 + 2x + 1$')
ax.text(x, y, r'$x_1, x_2, \ldots, x_n$')

# Greek letters
ax.text(x, y, r'$\alpha, \beta, \gamma, \delta, \epsilon$')
```

### 使用 Full LaTeX

```python
# Enable full LaTeX rendering (requires LaTeX installation)
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

ax.set_title(r'\textbf{Bold Title}')
ax.set_xlabel(r'Time $t$ (s)')
```

## 脊柱和网格

### 脊柱定制

```python
# Hide specific spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Move spine position
ax.spines['left'].set_position(('outward', 10))
ax.spines['bottom'].set_position(('data', 0))

# Change spine color and width
ax.spines['left'].set_color('red')
ax.spines['bottom'].set_linewidth(2)
```

### 网格定制

```python
# Basic grid
ax.grid(True)

# Customized grid
ax.grid(True, which='major', linestyle='--', linewidth=0.8, alpha=0.3)
ax.grid(True, which='minor', linestyle=':', linewidth=0.5, alpha=0.2)

# Grid for specific axis
ax.grid(True, axis='x')  # Only vertical lines
ax.grid(True, axis='y')  # Only horizontal lines

# Grid behind or in front of data
ax.set_axisbelow(True)  # Grid behind data
```

## 图例定制

### 图例定位

```python
# Location strings
ax.legend(loc='best')  # Automatic best position
ax.legend(loc='upper right')
ax.legend(loc='upper left')
ax.legend(loc='lower right')
ax.legend(loc='lower left')
ax.legend(loc='center')
ax.legend(loc='upper center')
ax.legend(loc='lower center')
ax.legend(loc='center left')
ax.legend(loc='center right')

# Precise positioning (bbox_to_anchor)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  # Outside plot area
ax.legend(bbox_to_anchor=(0.5, -0.15), loc='upper center', ncol=3)  # Below plot
```

### 传奇造型

```python
ax.legend(
    fontsize=12,
    frameon=True,           # Show frame
    framealpha=0.9,         # Frame transparency
    fancybox=True,          # Rounded corners
    shadow=True,            # Shadow effect
    ncol=2,                 # Number of columns
    title='Legend Title',   # Legend title
    title_fontsize=14,      # Title font size
    edgecolor='black',      # Frame edge color
    facecolor='white'       # Frame background color
)
```

### 自定义图例条目

```python
from matplotlib.lines import Line2D

# Create custom legend handles
custom_lines = [Line2D([0], [0], color='red', lw=2),
                Line2D([0], [0], color='blue', lw=2, linestyle='--'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10)]

ax.legend(custom_lines, ['Label 1', 'Label 2', 'Label 3'])
```

## 布局和间距

### 约束布局

```python
# Preferred method (automatic adjustment)
fig, axes = plt.subplots(2, 2, constrained_layout=True)
```

### 紧凑布局

```python
# Alternative method
fig, axes = plt.subplots(2, 2)
plt.tight_layout(pad=1.5, h_pad=2.0, w_pad=2.0)
```

### 手动调整

```python
# Fine-grained control
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1,
                    hspace=0.3, wspace=0.4)
```

## 专业的出版风格

出版物质量数据的示例配置：

```python
# Publication style configuration
plt.rcParams.update({
    # Figure
    'figure.figsize': (8, 6),
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,

    # Font
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 11,

    # Axes
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'axes.linewidth': 1.5,
    'axes.grid': False,
    'axes.spines.top': False,
    'axes.spines.right': False,

    # Lines
    'lines.linewidth': 2,
    'lines.markersize': 8,

    # Ticks
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',

    # Legend
    'legend.fontsize': 10,
    'legend.frameon': True,
    'legend.framealpha': 1.0,
    'legend.edgecolor': 'black'
})
```

## 黑暗主题

```python
# Dark background style
plt.style.use('dark_background')

# Or manual configuration
plt.rcParams.update({
    'figure.facecolor': '#1e1e1e',
    'axes.facecolor': '#1e1e1e',
    'axes.edgecolor': 'white',
    'axes.labelcolor': 'white',
    'text.color': 'white',
    'xtick.color': 'white',
    'ytick.color': 'white',
    'grid.color': 'gray',
    'legend.facecolor': '#1e1e1e',
    'legend.edgecolor': 'white'
})
```

## 颜色辅助功能

### 色盲友好调色板

```python
# Use colorblind-friendly colormaps
colorblind_friendly = ['viridis', 'plasma', 'cividis']

# Colorblind-friendly discrete colors
cb_colors = ['#0173B2', '#DE8F05', '#029E73', '#CC78BC',
             '#CA9161', '#949494', '#ECE133', '#56B4E9']

# Test with simulation tools or use these validated palettes
```

### 高对比度

```python
# Ensure sufficient contrast
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
```