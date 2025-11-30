<!-- 此文件由机器翻译自 operators.md -->

# Pymoo 遗传算子参考

pymoo 中遗传算子的综合参考。

## 采样运算符

抽样算子在优化开始时初始化总体。

### 随机抽样
**目的：**生成随机初始解
**类型：**
- `FloatRandomSampling`：连续变量
- `BinaryRandomSampling`：二进制变量
- `IntegerRandomSampling`：整数变量
- `PermutationRandomSampling`：基于排列的问题

**用途：**
```python
from pymoo.operators.sampling.rnd import FloatRandomSampling
sampling = FloatRandomSampling()
```

### 拉丁超立方采样 (LHS)
**目的：** 空间填充初始种群
**优点：** 比随机搜索空间更好的覆盖范围
**类型：**
- `LHS`：标准拉丁超立方体

**用途：**
<<<代码块_1>>>

### 定制采样
通过 Population 对象或 NumPy 数组提供初始种群

## 选择运算符

选择算子选择父母进行繁殖。

### 比赛选择
**目的：**通过锦标赛选拔家长
**机制：** 随机选择k个人，选择最好的
**参数：**
- `pressure`：锦标赛规模（默认值：2）
- `func_comp`：比较函数

**用途：**
<<<代码块_2>>>

### 随机选择
**目的：** 统一随机选择亲本
**用例：** 基线或探索型算法

**用途：**
<<<代码块_3>>>

## 交叉算子

交叉算子重新组合父解决方案以创建后代。

### 对于连续变量

#### 模拟二元交叉 (SBX)
**目的：** 持续优化的主要交叉
**机制：** 模拟二进制编码变量的单点交叉
**参数：**
- `prob`：交叉概率（默认值：0.9）
- `eta`：分布索引（默认值：15）
  - 更高的eta → 后代更接近父母
  - 更低的eta→更多的探索

**用途：**
<<<代码块_4>>>

**字符串简写：** `"real_sbx"`

#### 差分进化交叉
**目的：** DE 特异性重组
**变体：**
- `DE/rand/1/bin`
- `DE/best/1/bin`
- `DE/current-to-best/1/bin`

**参数：**
- `CR`：交叉率
- `F`：缩放因子

### 对于二元变量

#### 单点交叉
**目的：** 在一点处剪切和交换
**用途：**
<<<代码块_5>>>

#### 两点交叉
**用途：** 两点之间的剪切和交换
**用途：**
<<<代码块_6>>>

#### K 点交叉
**用途：** 多个切点
**参数：**
- `n_points`：交叉点的数量

#### 统一交叉
**目的：** 每个基因独立于亲本
**参数：**
- `prob`：每基因交换概率（默认值：0.5）

**用途：**
```python
from pymoo.operators.crossover.ux import UniformCrossover
crossover = UniformCrossover(prob=0.5)
```

#### 半均匀交叉 (HUX)
**目的：** 交换恰好一半的不同基因
**好处：** 保持遗传多样性

### 对于排列

#### 订单交叉 (OX)
**目的：** 维护父母的相对秩序
**用例：** 旅行推销员、日程安排问题

**用途：**
```python
from pymoo.operators.crossover.ox import OrderCrossover
crossover = OrderCrossover()
```

#### 边缘重组交叉 (ERX)
**目的：**保留父母的边缘信息
**用例：** 边缘连接很重要的路由问题

#### 部分映射交叉 (PMX)
**目的：** 交换段，同时保持排列有效性

## 变异算子

变异算子引入变异来保持多样性。

### 对于连续变量

#### 多项式变异 (PM)
**目的：** 用于持续优化的初级突变
**机制：** 多项式概率分布
**参数：**
- `prob`：每个变量的突变概率
- `eta`：分布索引（默认值：20）
  - 更高的eta→更小的扰动
  - 较低的eta→较大的扰动

**用途：**
```python
from pymoo.operators.mutation.pm import PM
mutation = PM(prob=None, eta=20)  # prob=None means 1/n_var
```

**字符串简写：** `"real_pm"`

**概率指南：**
- `None` 或 `1/n_var`：标准推荐
- 更高，更多探索
- 降低更多剥削

### 对于二元变量

#### 位翻转突变
**目的：** 以指定概率翻转位
**参数：**
- `prob`：每位翻转概率

**用途：**
```python
from pymoo.operators.mutation.bitflip import BitflipMutation
mutation = BitflipMutation(prob=0.05)
```

### 对于整数变量

#### 整数多项式变异
**目的：** PM 适用于整数
**确保：** 突变后的有效整数值

### 对于排列

#### 反转突变
**目的：** 反转一段排列
**用例：** 维护一些订单结构

**用途：**
```python
from pymoo.operators.mutation.inversion import InversionMutation
mutation = InversionMutation()
```
#### 扰乱突变
**目的：** 随机打乱一个片段

### 自定义突变
通过扩展 `Mutation` 类定义自定义突变

## 维修操作员

维修操作员修复约束违规或确保解决方案的可行性。

### 舍入修复
**用途：**四舍五入到最接近的有效值
**用例：** 具有边界约束的整数/离散变量

### 反弹修复
**目的：** 将越界值反映回可行区域
**用例：** 盒子约束的连续问题

### 投影修复
**目的：** 将不可行的解决方案投影到可行区域上
**用例：** 线性约束

### 定制维修
**目的：** 特定领域的约束处理
**实现：** 扩展 `Repair` 类

**示例：**
```python
from pymoo.core.repair import Repair

class MyRepair(Repair):
    def _do(self, problem, X, **kwargs):
        # Modify X to satisfy constraints
        # Return repaired X
        return X
```

## 操作员配置指南

### 参数调整

**交叉概率：**
- 高（0.8-0.95）：大多数问题的标准
- 较低：更强调突变

**突变概率：**
- `1/n_var`：标准推荐
- 更高：探索更多，收敛速度更慢
- 较低：收敛速度更快，有过早收敛的风险

**分布指数（eta）：**
- 交叉 eta (15-30)：本地搜索更高
- 突变 eta (20-50)：更高的开发利用

### 针对特定问题的选择

**持续出现的问题：**
- 交叉：SBX
- 变异：多项式变异
- 选择：锦标赛

**二元问题：**
- 交叉：两点或均匀
- 突变：Bitflip
- 选择：锦标赛

**排列问题：**
- 交叉：订单交叉（OX）
- 突变：倒置或扰乱
- 选择：锦标赛

**混合变量问题：**
- 每个变量类型使用适当的运算符
- 确保操作员兼容性

### 基于字符串的配置

Pymoo 支持方便的基于字符串的运算符规范：

```python
from pymoo.algorithms.soo.nonconvex.ga import GA

algorithm = GA(
    pop_size=100,
    sampling="real_random",
    crossover="real_sbx",
    mutation="real_pm"
)
```

**可用字符串：**
- 采样：`"real_random"`、`"real_lhs"`、`"bin_random"`、`"perm_random"`
- 交叉：`"real_sbx"`、`"real_de"`、`"int_sbx"`、`"bin_ux"`、`"bin_hux"`
- 突变：`"real_pm"`、`"int_pm"`、`"bin_bitflip"`、`"perm_inv"`

## 运算符组合示例

### 标准连续遗传算法：
```python
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.selection.tournament import TournamentSelection

sampling = FloatRandomSampling()
crossover = SBX(prob=0.9, eta=15)
mutation = PM(eta=20)
selection = TournamentSelection()
```

### 二进制遗传算法：
```python
from pymoo.operators.sampling.rnd import BinaryRandomSampling
from pymoo.operators.crossover.pntx import TwoPointCrossover
from pymoo.operators.mutation.bitflip import BitflipMutation

sampling = BinaryRandomSampling()
crossover = TwoPointCrossover()
mutation = BitflipMutation(prob=0.05)
```

### 排列 GA (TSP)：
```python
from pymoo.operators.sampling.rnd import PermutationRandomSampling
from pymoo.operators.crossover.ox import OrderCrossover
from pymoo.operators.mutation.inversion import InversionMutation

sampling = PermutationRandomSampling()
crossover = OrderCrossover()
mutation = InversionMutation()
```