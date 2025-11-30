<!-- 此文件由机器翻译自 function_reference.md -->

# Seaborn 函数参考

本文档按类别提供了所有主要 Seaborn 功能的综合参考。

## 关系图

### 散点图()

**目的：** 创建散点图，其中的点代表各个观察值。

**关键参数：**
- `data` - DataFrame、数组或数组字典
- `x, y` - x 和 y 轴的变量
- `hue` - 用于颜色编码的分组变量
- `size` - 用于大小编码的分组变量
- `style` - 标记样式的分组变量
- `palette` - 调色板名称或列表
- `hue_order` - 分类色调级别的顺序
- `hue_norm` - 数字色调标准化（元组或标准化对象）
- `sizes` - 大小编码的大小范围（元组或字典）
- `size_order` - 分类大小级别的顺序
- `size_norm` - 数字大小标准化
- `markers` - 标记样式（字符串、列表或字典）
- `style_order` - 分类样式级别的顺序
- `legend` - 如何绘制图例：“auto”、“brief”、“full”或 False
- `ax` - 要绘制的 Matplotlib 轴

**示例：**
```python
sns.scatterplot(data=df, x='height', y='weight',
                hue='gender', size='age', style='smoker',
                palette='Set2', sizes=(20, 200))
```

### 线图()

**目的：** 绘制具有自动聚合和重复测量置信区间的线图。

**关键参数：**
- `data` - DataFrame、数组或数组字典
- `x, y` - x 和 y 轴变量
- `hue` - 用于颜色编码的分组变量
- `size` - 线宽分组变量
- `style` - 线条样式（破折号）的分组变量
- `units` - 抽样单位的分组变量（单位内无聚合）
- `estimator` - 聚合观察值的函数（默认值：平均值）
- `errorbar` - 误差线方法：“sd”、“se”、“pi”、(“ci”、level)、(“pi”、level) 或 None
- `n_boot` - CI 计算的引导迭代次数
- `seed` - 用于可重复引导的随机种子
- `sort` - 绘图前对数据进行排序
- `err_style` - 用于错误表示的“band”或“bars”
- `err_kws` - 用于错误表示的附加参数
- `markers` - 用于强调数据点的标记样式
- `dashes` - 线条的虚线样式
- `legend` - 如何绘制图例
- `ax` - 要绘制的 Matplotlib 轴

**示例：**
<<<代码块_1>>>

### 相关图()

**用途：** 用于在 FacetGrid 上绘制关系图（散点图或线图）的图形级界面。

**关键参数：**
来自 `scatterplot()` 和 `lineplot()` 的所有参数，加上：
- `kind` - “分散”或“线条”
- `col` - 列面的分类变量
- `row` - 行面的分类变量
- `col_wrap` - 在这么多列之后换行
- `col_order` - 列构面级别的顺序
- `row_order` - 行构面级别的顺序
- `height` - 每个面的高度（以英寸为单位）
- `aspect` - 纵横比（宽度 = 高度 * 纵横比）
- `facet_kws` - FacetGrid 的附加参数

**示例：**
<<<代码块_2>>>

## 分布图

### 直方图()

**目的：** 通过灵活的分箱绘制单变量或双变量直方图。

**关键参数：**
- `data` - DataFrame、数组或字典
- `x, y` - 变量（y 对于二元变量是可选的）
- `hue` - 分组变量
- `weights` - 用于加权观察的变量
- `stat` - 聚合统计：“计数”、“频率”、“概率”、“百分比”、“密度”
- `bins` - bin 数量、bin 边缘或方法（“auto”、“fd”、“doane”、“scott”、“stone”、“rice”、“sturges”、“sqrt”）
- `binwidth` - bin 的宽度（覆盖 bin）
- `binrange` - 分箱范围（元组）
- `discrete` - 将 x 视为离散（将条形图置于值的中心）
- `cumulative` - 计算累积分布
- `common_bins` - 对所有色调级别使用相同的容器
- `common_norm` - 跨色调级别标准化
- `multiple` - 如何处理色调：“layer”、“dodge”、“stack”、“fill”
- `element` - 视觉元素：“bars”、“step”、“poly”
- `fill` - 填充条/元素
- `shrink` - 比例尺宽度（对于 multiple="dodge"）
- `kde` - 覆盖 KDE 估计
- `kde_kws` - KDE 参数
- `line_kws` - 步骤/多边形元素的参数
- `thresh` - bin 的最小计数阈值
- `pthresh` - 最小概率阈值
- `pmax` - 颜色缩放的最大概率
- `log_scale` - 轴的对数刻度（布尔值或基数）
- `legend` - 是否显示图例
- `ax` - Matplotlib 轴

**示例：**
<<<代码块_3>>>

### kdeplot()

**目的：** 绘制单变量或双变量核密度估计值。

**关键参数：**
- `data` - DataFrame、数组或字典
- `x, y` - 变量（y 对于二元变量是可选的）
- `hue` - 分组变量
- `weights` - 用于加权观察的变量
- `palette` - 调色板
- `hue_order` - 色调级别的顺序
- `hue_norm` - 数字色调标准化
- `multiple` - 如何处理色调：“layer”、“stack”、“fill”
- `common_norm` - 跨色调级别标准化
- `common_grid` - 对所有色调级别使用相同的网格
- `cumulative` - 计算累积分布
- `bw_method` - 带宽方法：“scott”、“silverman”或标量
- `bw_adjust` - 带宽乘数（越高=越平滑）
- `log_scale` - 轴的对数刻度
- `levels` - 轮廓级别的数字或值（双变量）
- `thresh` - 轮廓的最小密度阈值
- `gridsize` - 网格分辨率
- `cut` - 超出数据极限的扩展（以带宽为单位）
- `clip` - 曲线的数据范围（元组）
- `fill` - 填充曲线/轮廓下的区域
- `legend` - 是否显示图例
- `ax` - Matplotlib 轴

**示例：**
<<<代码块_4>>>

### ecdfplot()

**目的：** 绘制经验累积分布函数。

**关键参数：**
- `data` - DataFrame、数组或字典
- `x, y` - 变量（指定一个）
- `hue` - 分组变量
- `weights` - 用于加权观察的变量
- `stat` - “比例”或“计数”
- `complementary` - 绘制互补 CDF (1 - ECDF)
- `palette` - 调色板
- `hue_order` - 色调级别的顺序
- `hue_norm` - 数字色调标准化
- `log_scale` - 轴的对数刻度
- `legend` - 是否显示图例
- `ax` - Matplotlib 轴

**示例：**
<<<代码块_5>>>

### rugplot()

**目的：** 绘制刻度线，显示沿轴的各个观察结果。

**关键参数：**
- `data` - DataFrame、数组或字典
- `x, y` - 变量（指定一个）
- `hue` - 分组变量
- `height` - 刻度线高度（轴的比例）
- `expand_margins` - 为地毯添加边距空间
- `palette` - 调色板
- `hue_order` - 色调级别的顺序
- `hue_norm` - 数字色调标准化
- `legend` - 是否显示图例
- `ax` - Matplotlib 轴

**示例：**
<<<代码块_6>>>

### 分布图()

**用途：** 用于在 FacetGrid 上绘制分布图的图形级界面。

**关键参数：**
来自 `histplot()`、`kdeplot()` 和 `ecdfplot()` 的所有参数，加上：
- `kind` - “hist”、“kde”、“ecdf”
- `rug` - 在边缘轴上添加地毯图
- `rug_kws` - 地毯图参数
- `col` - 列面的分类变量
- `row` - 行面的分类变量
- `col_wrap` - 换行列
- `col_order` - 列面的顺序
- `row_order` - 行面的顺序
- `height` - 每个面的高度
- `aspect` - 纵横比
- `facet_kws` - FacetGrid 的附加参数

**示例：**
```python
sns.displot(data=df, x='measurement', hue='treatment',
            col='timepoint', kind='kde', fill=True,
            height=3, aspect=1.5, rug=True)
```

### 联合图()

**目的：** 使用边际单变量图绘制双变量图。

**关键参数：**
- `data` - 数据帧
- `x, y` - x 和 y 轴变量
- `hue` - 分组变量
- `kind` - “分散”、“kde”、“hist”、“十六进制”、“reg”、“resid”
- `height` - 图形大小（正方形）
- `ratio` - 关节与边缘轴的比率
- `space` - 关节轴和边缘轴之间的空间
- `dropna` - 删除缺失值
- `xlim, ylim` - 轴限制（元组）
- `marginal_ticks` - 在边缘轴上显示刻度
- `joint_kws` - 联合图参数
- `marginal_kws` - 边缘图的参数
- `hue_order` - 色调级别的顺序
- `palette` - 调色板

**示例：**
```python
sns.jointplot(data=df, x='var1', y='var2', hue='group',
              kind='scatter', height=6, ratio=4,
              joint_kws={'alpha': 0.5})
```

### 配对图()

**目的：** 绘制数据集中的成对关系。

**关键参数：**
- `data` - 数据帧
- `hue` - 用于颜色编码的分组变量
- `hue_order` - 色调级别的顺序
- `palette` - 调色板
- `vars` - 要绘制的变量（默认值：全数字）
- `x_vars, y_vars` - x 和 y 轴的变量（非方形网格）
- `kind` - “分散”、“kde”、“hist”、“reg”
- `diag_kind` - “自动”、“hist”、“kde”、无
- `markers` - 标记样式
- `height` - 每个面的高度
- `aspect` - 纵横比
- `corner` - 仅绘制下三角形
- `dropna` - 删除缺失值
- `plot_kws` - 非对角线图的参数
- `diag_kws` - 对角线图的参数
- `grid_kws` - PairGrid 的参数

**示例：**
```python
sns.pairplot(data=df, hue='species', palette='Set2',
             vars=['sepal_length', 'sepal_width', 'petal_length'],
             corner=True, height=2.5)
```

## 分类图

### 条形图()

**目的：** 绘制带有抖动点的分类散点图。

**关键参数：**
- `data` - DataFrame、数组或字典
- `x, y` - 变量（一个校准类别，一个连续变量）
- `hue` - 分组变量
- `order` - 分类级别的顺序
- `hue_order` - 色调级别的顺序
- `jitter` - 抖动量：True、float 或 False
- `dodge` - 并排单独的色调级别
- `orient` - “v”或“h”（通常推断）
- `color` - 所有元素均为单色
- `palette` - 调色板
- `size` - 标记大小
- `edgecolor` - 标记边缘颜色
- `linewidth` - 标记边缘宽度
- `native_scale` - 使用数字刻度作为分类轴
- `formatter` - 分类轴的格式化程序
- `legend` - 是否显示图例
- `ax` - Matplotlib 轴

**示例：**
```python
sns.stripplot(data=df, x='day', y='total_bill',
              hue='sex', dodge=True, jitter=0.2)
```

### 群图()

**目的：** 绘制具有不重叠点的分类散点图。

**关键参数：**
与 `stripplot()` 相同，但以下情况除外：
- 没有 `jitter` 参数
- `size` - 标记大小（对于避免重叠很重要）
- `warn_thresh` - 警告点过多的阈值（默认值：0.05）

**注意：** 对于大型数据集，计算量很大。对 >1000 个点使用条形图。

**示例：**
```python
sns.swarmplot(data=df, x='day', y='total_bill',
              hue='time', dodge=True, size=5)
```

### 箱线图()

**目的：** 绘制显示四分位数和异常值的箱线图。

**关键参数：**
- `data` - DataFrame、数组或字典
- `x, y` - 变量（一个校准类别，一个连续变量）
- `hue` - 分组变量
- `order` - 分类级别的顺序
- `hue_order` - 色调级别的顺序
- `orient` - “v”或“h”
- `color` - 盒子单色
- `palette` - 调色板
- `saturation` - 颜色饱和度强度
- `width` - 框的宽度
- `dodge` - 并排单独的色调级别
- `fliersize` - 离群值标记的大小
- `linewidth` - 方框线宽
- `whis` - 晶须的 IQR 乘数（默认值：1.5）
- `notch` - 绘制缺口框
- `showcaps` - 显示胡须帽
- `showmeans` - 显示平均值
- `meanprops` - 平均标记的属性
- `boxprops` - 框的属性
- `whiskerprops` - 晶须属性
- `capprops` - 大写字母的属性
- `flierprops` - 异常值的属性
- `medianprops` - 中线的属性
- `native_scale` - 使用数字刻度
- `formatter` - 分类轴的格式化程序
- `legend` - 是否显示图例
- `ax` - Matplotlib 轴

**示例：**
```python
sns.boxplot(data=df, x='day', y='total_bill',
            hue='smoker', palette='Set3',
            showmeans=True, notch=True)
```

### 小提琴图()

**目的：** 结合箱线图和 KDE 绘制小提琴图。

**关键参数：**
与 `boxplot()` 相同，加上：
- `bw_method` - KDE 带宽方法
- `bw_adjust` - KDE 带宽倍增器
- `cut` - KDE 超越极限的扩展
- `density_norm` - “面积”、“计数”、“宽度”
- `inner` - “盒子”、“四分位数”、“点”、“棒”、无
- `split` - 拆分小提琴以进行色调比较
- `scale` - 缩放方法：“面积”、“计数”、“宽度”
- `scale_hue` - 跨色调级别缩放
- `gridsize` - KDE 网格分辨率

**示例：**
```python
sns.violinplot(data=df, x='day', y='total_bill',
               hue='sex', split=True, inner='quartile',
               palette='muted')
```

### 箱形图()

**目的：** 为显示更多分位数的较大数据集绘制增强的箱线图。

**关键参数：**
与 `boxplot()` 相同，加上：
- `k_depth` - “tukey”、“比例”、“可信”、“完整”或 int
- `outlier_prop` - 异常值数据的比例
- `trust_alpha` - Alpha 值得信赖的深度
- `showfliers` - 显示异常点

**示例：**
```python
sns.boxenplot(data=df, x='day', y='total_bill',
              hue='time', palette='Set2')
```

### 条形图()

**目的：** 绘制带有误差线的条形图，显示统计估计值。

**关键参数：**
- `data` - DataFrame、数组或字典
- `x, y` - 变量（一个校准类别，一个连续变量）
- `hue` - 分组变量
- `order` - 分类级别的顺序
- `hue_order` - 色调级别的顺序
- `estimator` - 聚合函数（默认值：均值）
- `errorbar` - 错误表示：“sd”、“se”、“pi”、(“ci”、level)、(“pi”、level) 或 None
- `n_boot` - 引导迭代
- `seed` - 随机种子
- `units` - 采样单位标识符
- `weights` - 观察权重
- `orient` - “v”或“h”
- `color` - 单条颜色
- `palette` - 调色板
- `saturation` - 颜色饱和度
- `width` - 条形宽度
- `dodge` - 并排单独的色调级别
- `errcolor` - 错误栏颜色
- `errwidth` - 误差线线宽
- `capsize` - 错误栏上限宽度
- `native_scale` - 使用数字刻度
- `formatter` - 分类轴的格式化程序
- `legend` - 是否显示图例
- `ax` - Matplotlib 轴

**示例：**
```python
sns.barplot(data=df, x='day', y='total_bill',
            hue='sex', estimator='median',
            errorbar=('ci', 95), capsize=0.1)
```

### 计数图()

**目的：** 显示每个分类箱中的观察计数。

**关键参数：**
与 `barplot()` 相同，但是：
- 仅指定 x 或 y 之一（分类变量）
- 没有估计器或误差条（显示计数）
- `stat` - “计数”或“百分比”

**示例：**
```python
sns.countplot(data=df, x='day', hue='time',
              palette='pastel', dodge=True)
```

### 点图()

**目的：** 用连接线显示点估计值和置信区间。

**关键参数：**
与 `barplot()` 相同，加上：
- `markers` - 标记样式
- `linestyles` - 线条样式
- `scale` - 标记比例
- `join` - 用线连接点
- `capsize` - 错误栏上限宽度

**示例：**
```python
sns.pointplot(data=df, x='time', y='total_bill',
              hue='sex', markers=['o', 's'],
              linestyles=['-', '--'], capsize=0.1)
```

###猫图()

**用途：** FacetGrid 上分类图的图形级界面。

**关键参数：**
分类图中的所有参数，加上：
- `kind` - “条”、“群”、“盒子”、“小提琴”、“盒子”、“条”、“点”、“计数”
- `col` - 列面的分类变量
- `row` - 行面的分类变量
- `col_wrap` - 换行列
- `col_order` - 列面的顺序
- `row_order` - 行面的顺序
- `height` - 每个面的高度
- `aspect` - 纵横比
- `sharex, sharey` - 跨面共享轴
- `legend` - 是否显示图例
- `legend_out` - 将图例放在图之外
- `facet_kws` - 其他 FacetGrid 参数

**示例：**
```python
sns.catplot(data=df, x='day', y='total_bill',
            hue='smoker', col='time',
            kind='violin', split=True,
            height=4, aspect=0.8)
```

## 回归图

### regplot()

**目的：** 绘制数据和线性回归拟合。

**关键参数：**
- `data` - 数据帧
- `x, y` - 变量或数据向量
- `x_estimator` - 将估计器应用于 x bin
- `x_bins` - 用于估计器的 Bin x
- `x_ci` - 用于分箱估计的 CI
- `scatter` - 显示散点
- `fit_reg` - 绘制回归线
- `ci` - 用于回归估计的 CI（int 或 None）
- `n_boot` - CI 的引导迭代
- `units` - 采样单位标识符
- `seed` - 随机种子
- `order` - 多项式回归阶数
- `logistic` - 拟合逻辑回归
- `lowess` - 使 Lowess 更平滑
- `robust` - 拟合稳健回归
- `logx` - 对数变换 x
- `x_partial, y_partial` - 部分回归（回归变量）
- `truncate` - 将回归线限制在数据范围内
- `dropna` - 删除缺失值
- `x_jitter, y_jitter` - 向数据添加抖动
- `label` - 图例标签
- `color` - 所有元素的颜色
- `marker` - 标记样式
- `scatter_kws` - 分散参数
- `line_kws` - 回归线参数
- `ax` - Matplotlib 轴

**示例：**
```python
sns.regplot(data=df, x='total_bill', y='tip',
            order=2, robust=True, ci=95,
            scatter_kws={'alpha': 0.5})
```

### lmplot()

**用途：** 用于在 FacetGrid 上绘制回归图的图形级界面。

**关键参数：**
来自 `regplot()` 的所有参数，加上：
- `hue` - 分组变量
- `col` - 列面
- `row` - 行面
- `palette` - 调色板
- `col_wrap` - 换行列
- `height` - 小平面高度
- `aspect` - 纵横比
- `markers` - 标记样式
- `sharex, sharey` - 共享轴
- `hue_order` - 色调级别的顺序
- `col_order` - 列面的顺序
- `row_order` - 行面的顺序
- `legend` - 是否显示图例
- `legend_out` - 将图例放在外面
- `facet_kws` - FacetGrid 参数

**示例：**
```python
sns.lmplot(data=df, x='total_bill', y='tip',
           hue='smoker', col='time', row='sex',
           height=3, aspect=1.2, ci=None)
```

### 剩余图()

**目的：** 绘制回归的残差。

**关键参数：**
与 `regplot()` 相同，但是：
- 始终绘制残差（y - 预测）与 x 的关系图
- 在 y=0 处添加水平线
- `lowess` - 将 Lowess 平滑器拟合到残差

**示例：**
```python
sns.residplot(data=df, x='x', y='y', lowess=True,
              scatter_kws={'alpha': 0.5})
```

## 矩阵图

### 热图()

**目的：** 将矩形数据绘制为颜色编码矩阵。

**关键参数：**
- `data` - 二维数组数据
- `vmin, vmax` - 颜色图的锚点值
- `cmap` - 颜色图名称或对象
- `center` - 颜色图中心的值
- `robust` - 对颜色图范围使用稳健的分位数
- `annot` - 注释单元格：True、False 或数组
- `fmt` - 注释的格式字符串（例如“.2f”）
- `annot_kws` - 注释参数
- `linewidths` - 单元格边框的宽度
- `linecolor` - 单元格边框的颜色
- `cbar` - 绘制颜色条
- `cbar_kws` - 颜色条参数
- `cbar_ax` - 颜色条轴
- `square` - 强制方形单元格
- `xticklabels, yticklabels` - 勾选标签（True、False、int 或列表）
- `mask` - 用于屏蔽单元格的布尔数组
- `ax` - Matplotlib 轴

**示例：**
```python
# Correlation matrix
corr = df.corr()
mask = np.triu(np.ones_like(corr, dtype=bool))
sns.heatmap(corr, mask=mask, annot=True, fmt='.2f',
            cmap='coolwarm', center=0, square=True,
            linewidths=1, cbar_kws={'shrink': 0.8})
```

### 聚类图()

**目的：** 绘制分层聚类热图。

**关键参数：**
来自 `heatmap()` 的所有参数，加上：
- `pivot_kws` - 用于旋转的参数（如果需要）
- `method` - 联动方式：“单一”、“完整”、“平均”、“加权”、“质心”、“中值”、“病房”
- `metric` - 聚类的距离度量
- `standard_scale` - 标准化数据：0（行）、1（列）或无
- `z_score` - Z 分数标准化数据：0（行）、1（列）或无
- `row_cluster, col_cluster` - 集群行/列
- `row_linkage, col_linkage` - 预先计算的链接矩阵
- `row_colors, col_colors` - 附加颜色注释
- `dendrogram_ratio` - 树状图与热图的比率
- `colors_ratio` - 颜色注释与热图的比率
- `cbar_pos` - 颜色条位置（元组：x、y、宽度、高度）
- `tree_kws` - 树状图参数
- `figsize` - 图形尺寸

**示例：**
```python
sns.clustermap(data, method='average', metric='euclidean',
               z_score=0, cmap='viridis',
               row_colors=row_colors, col_colors=col_colors,
               figsize=(12, 12), dendrogram_ratio=0.1)
```

## 多图网格

### FacetGrid

**用途：** 用于绘制条件关系的多图网格。

**初始化：**
```python
g = sns.FacetGrid(data, row=None, col=None, hue=None,
                  col_wrap=None, sharex=True, sharey=True,
                  height=3, aspect=1, palette=None,
                  row_order=None, col_order=None, hue_order=None,
                  hue_kws=None, dropna=False, legend_out=True,
                  despine=True, margin_titles=False,
                  xlim=None, ylim=None, subplot_kws=None,
                  gridspec_kws=None)
```

**方法：**
- `map(func, *args, **kwargs)` - 将函数应用于每个方面
- `map_dataframe(func, *args, **kwargs)` - 使用完整的 DataFrame 应用函数
- `set_axis_labels(x_var, y_var)` - 设置轴标签
- `set_titles(template, **kwargs)` - 设置子图标题
- `set(kwargs)` - 在所有轴上设置属性
- `add_legend(legend_data, title, label_order, **kwargs)` - 添加图例
- `savefig(*args, **kwargs)` - 保存图形

**示例：**
```python
g = sns.FacetGrid(df, col='time', row='sex', hue='smoker',
                  height=3, aspect=1.5, margin_titles=True)
g.map(sns.scatterplot, 'total_bill', 'tip', alpha=0.7)
g.add_legend()
g.set_axis_labels('Total Bill ($)', 'Tip ($)')
g.set_titles('{col_name} | {row_name}')
```

### 配对网格

**用途：** 用于绘制数据集中的成对关系的网格。

**初始化：**
```python
g = sns.PairGrid(data, hue=None, vars=None,
                 x_vars=None, y_vars=None,
                 hue_order=None, palette=None,
                 hue_kws=None, corner=False,
                 diag_sharey=True, height=2.5,
                 aspect=1, layout_pad=0.5,
                 despine=True, dropna=False)
```

**方法：**
- `map(func, **kwargs)` - 将函数应用于所有子图
- `map_diag(func, **kwargs)` - 应用于对角线
- `map_offdiag(func, **kwargs)` - 适用于非对角线
- `map_upper(func, **kwargs)` - 应用于上三角形
- `map_lower(func, **kwargs)` - 应用于下三角
- `add_legend(legend_data, **kwargs)` - 添加图例
- `savefig(*args, **kwargs)` - 保存图形

**示例：**
```python
g = sns.PairGrid(df, hue='species', vars=['a', 'b', 'c', 'd'],
                 corner=True, height=2.5)
g.map_upper(sns.scatterplot, alpha=0.5)
g.map_lower(sns.kdeplot)
g.map_diag(sns.histplot, kde=True)
g.add_legend()
```

### 联合网格

**目的：** 双变量图与边际单变量图的网格。

**初始化：**
```python
g = sns.JointGrid(data=None, x=None, y=None, hue=None,
                  height=6, ratio=5, space=0.2,
                  dropna=False, xlim=None, ylim=None,
                  marginal_ticks=False, hue_order=None,
                  palette=None)
```

**方法：**
- `plot(joint_func, marginal_func, **kwargs)` - 绘制关节和边缘
- `plot_joint(func, **kwargs)` - 绘制联合分布
- `plot_marginals(func, **kwargs)` - 绘制边际分布
- `refline(x, y, **kwargs)` - 添加参考线
- `set_axis_labels(xlabel, ylabel, **kwargs)` - 设置轴标签
- `savefig(*args, **kwargs)` - 保存图形

**示例：**
```python
g = sns.JointGrid(data=df, x='x', y='y', hue='group',
                  height=6, ratio=5, space=0.2)
g.plot_joint(sns.scatterplot, alpha=0.5)
g.plot_marginals(sns.histplot, kde=True)
g.set_axis_labels('Variable X', 'Variable Y')
```