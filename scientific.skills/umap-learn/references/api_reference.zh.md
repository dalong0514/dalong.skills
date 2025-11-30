<!-- 此文件由机器翻译自 api_reference.md -->

# UMAP API 参考

## UMAP 类

`umap.UMAP(n_neighbors=15, n_components=2, metric='euclidean', n_epochs=None, learning_rate=1.0, init='spectral', min_dist=0.1, spread=1.0, low_memory=True, set_op_mix_ratio=1.0, local_connectivity=1.0, repulsion_strength=1.0, negative_sample_rate=5, transform_queue_size=4.0, a=None, b=None, random_state=None, metric_kwds=None, angular_rp_forest=False, target_n_neighbors=-1, target_metric='categorical', target_metric_kwds=None, target_weight=0.5, transform_seed=42, transform_mode='embedding', force_approximation_algorithm=False, verbose=False, unique=False, densmap=False, dens_lambda=2.0, dens_frac=0.3, dens_var_shift=0.1, output_dens=False, disconnection_distance=None, precomputed_knn=(None, None, None))`

找到近似数据底层流形的低维嵌入。

### 核心参数

#### n_neighbors（整数，默认值：15）
用于流形近似的局部邻域的大小。较大的值会产生更多的流形全局视图，而较小的值会保留更多的局部结构。一般在2到100之间。

**调整指导：**
- 对于非常局部的结构使用 2-5
- 使用 10-20 平衡局部/全局结构（典型）
- 使用 50-200 强调全局结构

#### n_components（整数，默认值：2）
嵌入空间的尺寸。与 t-SNE 不同，UMAP 随着嵌入维度的增加而很好地扩展。

**共同价值观：**
- 2-3：可视化
- 5-10：聚类预处理
- 10-100：下游 ML 的特征工程

#### 公制（str 或可调用，默认值：'euclidean'）
要使用的距离度量。接受：
- 来自 scipy.spatial.distance 的任何度量
- sklearn.metrics 中的任何指标
- 自定义可调用距离函数（必须使用 Numba 编译）

**常用指标：**
- `'euclidean'`：标准欧氏距离（默认）
- `'manhattan'`：L1距离
- `'cosine'`：余弦距离（适用于文本/文档向量）
- `'correlation'`：相关距离
- `'hamming'`：汉明距离（对于二进制数据）
- `'jaccard'`：杰卡德距离（对于二进制/集合数据）
- `'dice'`：骰子距离
- `'canberra'`：堪培拉距离
- `'braycurtis'`：布雷-柯蒂斯距离
- `'chebyshev'`：切比雪夫距离
- `'minkowski'`：闵可夫斯基距离（用 metric_kwds 指定 p）
- `'precomputed'`：使用预先计算的距离矩阵

#### min_dist（浮点数，默认值：0.1）
嵌入点之间的有效最小距离。控制点堆积在一起的紧密程度。较小的值会导致更密集的嵌入。

**调整指导：**
- 对集群应用程序使用0.0
- 使用 0.1-0.3 进行可视化（平衡）
- 使用 0.5-0.99 保存松散结构

#### 点差（浮点数，默认值：1.0）
嵌入点的有效规模。与 `min_dist` 结合使用来控制聚集嵌入与分散嵌入。确定簇在嵌入空间中的分布情况。

### 训练参数

#### n_epochs（整数，默认值：无）
训练纪元数。如果无，则根据数据集大小（通常为 200-500 epoch）自动确定。

**手动调整：**
- 较小的数据集可能需要 500+ epoch
- 更大的数据集可能会在 200 个 epoch 内收敛
- 更多纪元 = 更好的优化，但训练速度更慢

####学习率（浮点数，默认值：1.0）
SGD 优化器的初始学习率。较高的值会导致更快的收敛，但可能会超出最佳解决方案。

#### init（str 或 np.ndarray，默认值：'spectral'）
嵌入的初始化方法：
- `'spectral'`：使用频谱嵌入（默认，通常是最好的）
- `'random'`：随机初始化
- `'pca'`：使用 PCA 初始化
- numpy 数组：自定义初始化（形状：(n_samples, n_components)）

### 高级结构参数

#### local_connectivity（整数，默认值：1.0）
假设本地连接的最近邻居的数量。值越高，连接的流形越多。

#### set_op_mix_ratio（浮点数，默认值：1.0）
构造模糊集并集时并集和交集之间的插值。值 1.0 使用纯并集，0.0 使用纯交集。

#### 斥力强度（浮点数，默认值：1.0）
低维嵌入优化中应用于负样本的加权。值越高，嵌入点的距离就越远。

#### 负样本率（整数，默认值：5）
每个正样本选择的负样本数。较高的值会导致点之间更大的排斥和更多的分散嵌入，但会增加计算成本。

### 监督学习参数

#### target_n_neighbors（整数，默认值：-1）
构造目标单纯集时要使用的最近邻数。如果为 -1，则使用 n_neighbors 值。

#### target_metric（str，默认值：“分类”）
目标值（标签）的距离度量：
- `'categorical'`：用于分类任务
- 回归任务的任何其他指标

#### target_weight（浮点数，默认值：0.5）
应用于目标信息与数据结构的权重。范围 0.0 至 1.0：
- 0.0：纯无监督嵌入（忽略标签）
- 0.5：平衡（默认）
- 1.0：纯监督嵌入（仅考虑标签）
### 变换参数

####transform_queue_size（浮点数，默认值：4.0）
用于变换操作的最近邻搜索队列的大小。较大的值可以提高转换精度，但会增加内存使用量和计算时间。

####transform_seed（整数，默认值：42）
用于变换操作的随机种子。确保转换结果的可重复性。

####transform_mode（str，默认值：'嵌入'）
转换新数据的方法：
- `'embedding'`：标准方法（默认）
- `'graph'`：使用最近邻图

### 性能参数

#### low_memory（布尔值，默认值：True）
是否使用内存高效的实现。仅当内存不受限制并且您想要更快的性能时才设置为 False。

#### 详细（布尔值，默认值：False）
是否在拟合过程中打印进度消息。

#### 唯一（布尔值，默认值：False）
是否仅考虑唯一数据点。如果您知道数据包含许多重复项，请设置为 True 以提高性能。

####force_approximation_algorithm（布尔值，默认值：False）
即使对于小型数据集，也强制使用近似最近邻搜索。可以提高大型数据集的性能。

#### angular_rp_forest（布尔值，默认值：False）
是否使用角度随机投影森林进行最近邻搜索。可以提高高维标准化数据的性能。

### DensMAP 参数

DensMAP 是保留局部密度信息的变体。

#### densmap（布尔值，默认值：False）
是否使用 DensMAP 算法而不是标准 UMAP。除了拓扑结构之外，还保留局部密度。

#### dens_lambda（浮点数，默认值：2.0）
DensMAP 优化中密度保留项的权重。较高的值强调密度保持。

#### dens_frac（浮点数，默认值：0.3）
用于 DensMAP 中密度估计的数据集的一部分。

#### dens_var_shift（浮点数，默认值：0.1）
DensMAP 中密度估计的正则化参数。

####output_dens（布尔值，默认值：False）
除了嵌入之外是否还输出局部密度估计。结果存储在 `rad_orig_` 和 `rad_emb_` 属性中。

### 其他参数

#### a（浮点数，默认值：无）
参数控制嵌入。如果没有，则根据 min_dist 和 spread 自动确定。

#### b（浮点数，默认值：无）
参数控制嵌入。如果没有，则根据 min_dist 和 spread 自动确定。

#### random_state（int、RandomState 实例或 None，默认值：None）
随机状态以实现可重复性。设置为整数以获得可重现的结果。

#### metric_kwds（字典，默认值：无）
距离度量的附加关键字参数。

####断开连接距离（浮点数，默认值：无）
考虑断开点的距离阈值。如果无，则使用图中的最大距离。

#### precompulated_knn（元组，默认值：（无，无，无））
预先计算的 k 最近邻为 (knn_indices, knn_dists, knn_search_index)。对于重用昂贵的计算很有用。

## 方法

### 适合（X，y=无）
将 UMAP 模型拟合到数据。

**参数：**
- `X`：类似数组，形状 (n_samples, n_features) - 训练数据
- `y`：类似数组，形状 (n_samples,)，可选 - 监督降维的目标值

**退货：**
- `self`：拟合的 UMAP 对象

**属性设置：**
- `embedding_`：训练数据的嵌入表示
- `graph_`：流形的模糊单纯形集逼近
- `_raw_data`：训练数据的副本
- `_small_data`：数据集是否被认为是小
- `_metric_kwds`：已处理的度量关键字参数
- `_n_neighbors`：实际使用的 n_neighbors
- `_initial_alpha`：初始学习率
- `_a`、`_b`：曲线参数

### fit_transform(X, y=无)
拟合模型并返回嵌入的表示。

**参数：**
- `X`：类似数组，形状 (n_samples, n_features) - 训练数据
- `y`：类似数组，形状 (n_samples,)，可选 - 监督降维的目标值

**退货：**
- `X_new`：数组，形状（n_samples，n_components） - 嵌入数据

### 变换(X)
将新数据转换到现有的嵌入空间中。

**参数：**
- `X`：类似数组，形状 (n_samples, n_features) - 要转换的新数据

**退货：**
- `X_new`：数组，形状（n_samples，n_components） - 新数据的嵌入表示

**重要说明：**
- 调用变换之前必须先拟合模型
- 转换质量取决于训练和测试分布之间的相似性
- 对于显着不同的数据分布，请考虑参数化 UMAP

### 逆变换(X)
将数据从嵌入空间变换回原始数据空间。

**参数：**
- `X`：类似数组，形状 (n_samples, n_components) - 嵌入数据点

**退货：**
- `X_new`：数组，形状（n_samples，n_features） - 原始空间中的重构数据

**重要说明：**
- 计算成本较高的操作
- 在训练嵌入的凸包之外效果不佳
- 重建质量因地区而异

### 更新(X)
使用新数据更新模型。允许增量拟合。

**参数：**
- `X`：类似数组，形状 (n_samples, n_features) - 要合并的新数据

**退货：**
- `self`：更新的 UMAP 对象

**注意：** 实验功能，可能不会保留批量训练的所有属性。

## 属性

### 嵌入_
array, shape (n_samples, n_components) - 训练数据的嵌入表示。

### 图_
scipy.sparse.csr_matrix - 流形的模糊单纯形集逼近的加权邻接矩阵。

### _raw_data
array - 原始训练数据的副本。

### _稀疏数据
bool - 训练数据是否稀疏。

### _小数据
bool - 数据集是否被视为小（对小数据集使用不同的算法）。

### _input_hash
str - 用于缓存目的的输入数据的哈希值。

### _knn_indices
array - 每个训练点的 k 最近邻的索引。

### _knn_dists
array - 每个训练点到 k 最近邻的距离。

### _rp_forest
list - 用于近似最近邻搜索的随机投影森林。

## 参数UMAP类

`umap.ParametricUMAP(encoder=None, decoder=None, parametric_reconstruction=False, autoencoder_loss=False, reconstruction_validation=None, dims=None, batch_size=None, n_training_epochs=1, loss_report_frequency=10, optimizer=None, keras_fit_kwargs={}, **kwargs)`

参数化 UMAP 使用神经网络来学习嵌入函数。

### 附加参数（UMAP 之外）

#### 编码器（tensorflow.keras.Model，默认值：无）
用于将数据编码为嵌入的 Keras 模型。如果无，则使用默认的 3 层架构，每层 100 个神经元。

#### 解码器（tensorflow.keras.Model，默认值：无）
用于将嵌入解码回数据空间的 Keras 模型。仅当 parametric_reconstruction=True 时使用。

#### parametric_reconstruction （布尔值，默认值：False）
是否使用参数化重建。需要解码器型号。

#### autoencoder_loss（布尔值，默认值：False）
优化中是否包含重建损失。需要解码器型号。

####重建_验证（元组，默认值：无）
用于监控训练期间重建损失的验证数据（X_val，y_val）。

#### 暗淡（元组，默认值：无）
编码器网络的输入维度。如果提供自定义编码器，则需要。

####batch_size（整数，默认值：无）
神经网络训练的批量大小。如果无，则自动确定。

#### n_training_epochs（整数，默认值：1）
神经网络的训练纪元数。更多纪元可以提高质量，但会增加训练时间。

#### loss_report_Frequency（整数，默认值：10）
培训期间报告损失的频率。

#### 优化器（tensorflow.keras.optimizers.Optimizer，默认值：无）
用于训练的 Keras 优化器。如果没有，则使用带有learning_rate参数的Adam。

#### keras_fit_kwargs（字典，默认值：{}）
传递给 Keras fit() 方法的其他关键字参数。

### 方法

与 UMAP 类相同，但 transform() 和 inverse_transform() 使用学习的神经网络来实现更快的推理。

## 实用函数

### umap.nearest_neighbors(X, n_neighbors, metric, metric_kwds={}, angle=False, random_state=None)
计算数据的 k 最近邻。

**返回：** (knn_indices, knn_dists, rp_forest)

### umap.fuzzy_simplicial_set(X, n_neighbors, random_state, metric, metric_kwds={}, knn_indices=None, knn_dists=None, Angular=False, set_op_mix_ratio=1.0, local_connectivity=1.0, apply_set_operations=True, verbose=False, return_dists=None)
构造数据的模糊单纯集表示。

**返回：** 模糊单纯集作为稀疏矩阵
### umap.simplicial_set_embedding（数据，图，n_components，initial_alpha，a，b，gamma，负采样率，n_epochs，init，random_state，metric，metric_kwds，densmap，densmap_kwds，output_dens，output_metric，output_metric_kwds，euclidean_output，parallel=False，verbose=False）
执行优化以找到低维嵌入。

**返回：** 嵌入数组

### umap.find_ab_params(spread, min_dist)
根据 spread 和 min_dist 拟合 UMAP 曲线的 a、b 参数。

**返回：** (a, b) 元组

## AlignedUMAP 类

`umap.AlignedUMAP(n_neighbors=15, n_components=2, metric='euclidean', alignment_regularisation=1e-2, alignment_window_size=3, **kwargs)`

用于对齐多个相关数据集的 UMAP 变体。

### 附加参数

####alignment_regularization（浮点数，默认值：1e-2）
数据集之间对齐正则化的强度。

####alignment_window_size（int，默认值：3）
要对齐的相邻数据集的数量。

### 方法

#### 适合（X）
将模型拟合到多个数据集。

**参数：**
- `X`：数组列表 - 要对齐的数据集列表

**退货：**
- `self`：拟合模型

### 属性

#### 嵌入_
数组列表 - 对齐嵌入的列表，每个输入数据集一个。

## 用法示例

### 所有常用参数的基本用法

```python
import umap

# Standard 2D visualization embedding
reducer = umap.UMAP(
    n_neighbors=15,          # Balance local/global structure
    n_components=2,          # Output dimensions
    metric='euclidean',      # Distance metric
    min_dist=0.1,           # Minimum distance between points
    spread=1.0,             # Scale of embedded points
    random_state=42,        # Reproducibility
    n_epochs=200,           # Training iterations (None = auto)
    learning_rate=1.0,      # SGD learning rate
    init='spectral',        # Initialization method
    low_memory=True,        # Memory-efficient mode
    verbose=True            # Print progress
)

embedding = reducer.fit_transform(data)
```

### 监督学习

<<<代码块_1>>>

### 聚类预处理

<<<代码块_2>>>

### 自定义距离度量

<<<代码块_3>>>

### 具有自定义架构的参数化 UMAP

<<<代码块_4>>>

### DensMAP 用于密度保持

<<<代码块_5>>>

### 时间序列的对齐 UMAP

<<<代码块_6>>>