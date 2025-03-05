# 文件说明
因为细胞具有批次效应，并且该效应比较明显，所以27号的分析将三个批次分开进行了预处理和聚类。

## 原始数据文件
原始数据文件为`Expr_Mat_PS+Mesoderm_Filter_ComBat_No.Rep_Cell.Name.Merge_2.txt`
测试数据文件为`test_data.txt`，只选取了原始数据文件中的前4个细胞用于代码测试，均来自第一批次。

## 结果文件
1. 三个批次预处理后的结果h5ad文件：`preprocessed_data_{batch_ID}.h5ad`
2. 三个批次预处理时产生的图表文件夹名：`preprocessed_figures_{batch_ID}`
3. 三个批次聚类后的结果h5ad文件：`clustered_data_{batch_ID}.h5ad`
4. 三个批次聚类时产生的图表文件夹名：`umap_clustering_figures_{batch_ID}`
5. 消除批次效应前，三个批次直接合并进行预处理的结果h5ad文件：`preprocessed_data_123.h5ad`
6. 消除批次效应前，三个批次直接合并进行预处理时产生的图表文件夹名：`preprocessed_figures_123`
7. 消除批次效应前，三个批次直接合并聚类后的结果h5ad文件：`clustered_data_123.h5ad`
8. 消除批次效应前，三个批次直接合并聚类时产生的图表文件夹名：`umap_clustering_figures_123`
9. `insight_figures`文件夹中的图片为`insight_n_genes_by_counts.py`结果文件。
10. `{}_noBatchEffect`后缀代表在`preprocess.py`消除批次效应后的所有结果。

## h5ad文件对obs和var的注释说明
### 对obs的注释
`batch`：表示细胞测序批次，值为 **1、2、3**，用于区分不同批次的实验数据。
`cell_type`：表示注释的细胞类别，用于标识每个细胞的类型。
`cell_time`：表示细胞发育时期，取值为 **E65、E675、D70、E725、E75**，用于描述细胞在发育过程中的时间点。
`cell_ID`：表示细胞的唯一标识符，用于区分不同的细胞。
`n_genes`：表示每个细胞中检测到的基因数量。
`n_genes_by_counts`：表示每个细胞中基于表达量（counts）检测到的基因数量。
`total_counts`：表示每个细胞的总表达量（counts），即所有基因表达量的总和。
`leiden`：表示基于 Leiden 聚类算法得到的细胞聚类标签，用于标识细胞所属的聚类群体。

### 对var的注释
`n_cells`：表示表达该基因的细胞数量。
`n_cells_by_counts`：表示基于表达量（counts）检测到该基因的细胞数量。
`mean_counts`：表示该基因在所有细胞中的平均表达量。
`pct_dropout_by_counts`：表示该基因在细胞中的丢失率（dropout rate），即未检测到该基因表达的细胞比例。
`total_counts`：表示该基因在所有细胞中的总表达量。
`highly_variable`：表示该基因是否被标记为高变基因（highly variable gene, HVG），True 表示是高变基因，False 表示不是。
`means`：表示该基因在标准化后的数据中的平均表达量。
`dispersions`：表示该基因的离散度（dispersion），用于识别高变基因。
`dispersions_norm`：表示该基因的归一化离散度（normalized dispersion），用于识别高变基因。
`mean`：表示该基因在标准化后的数据中的平均表达量（与 means 相同）。
`std`：表示该基因在标准化后的数据中的标准差。

## 脚本文件
### 脚本`preprocess.py`
#### 编写参考
按照scanpy的[标准数据处理流程处理](https://scanpy.readthedocs.io/en/1.11.x/tutorials/basics/clustering-2017.html#finding-marker-genes)

#### Tips
1. scanpy的工作流相对固定，**前后存在依赖**。请不要随意调整该脚本中函数的顺序。
2. 虽然聚类表征、后续的拟时间分析均采用了UMAP降维算法，但因为scanpy中执行umap降维工作流需要提前运行PCA算法，因此在该脚本中还出现了PCA函数，这不是多余的。
3. 数据`Expr_Mat_PS+Mesoderm_Filter_ComBat_No.Rep_Cell.Name.Merge_2.txt`经过老师预先处理，在清洗过程中全部细胞均通过质量控制。
4. 在`umap_clustering_figures_123`中可以观察到很明显的批次效应，三批的细胞`n_genes_by_counts`分布差别很大，分别分布在3k-4k，8-9.5k,10-11k的三个条带上。测序深度差别细微。
5. DeepSeek提到可以同时使用`sc.pp.regress_out`和`sc.pp.combat`来去除批次效应。并且提到“如果批次效应与某些技术性偏差（如测序深度）强相关，建议先去除这些技术性偏差，再处理批次效应”，而提供的数据表现出了明显的批次效应和测序深度存在较弱相关，两条建议脚本中均采用。
6. 在脚本中提取第一批次的细胞时，请写`adata_1 = adata[adata.obs.index.str.split('.').str[-2] == '1', :]`，而不能写成`adata_1 = adata[adata.obs.index.str.split('.')[-2] == 1, :]`。因为`.str[-2]`中的`str`是`pandas`提供的方法，代表对列表**逐元素操作，取每个元素的倒数第二个值**。

### 脚本`umap_clustering.py`
#### 编写参考
按照scanpy的[标准数据处理流程处理](https://scanpy.readthedocs.io/en/1.11.x/tutorials/basics/clustering-2017.html#finding-marker-genes)

#### 内容
包含对单细胞数据进行降维的代码，使用UMAP算法。但UMAP降维效果不好，算法没有学习到明显的流型。

### 脚本`show_batch_effect.py`
用于展示批次效应可能带来的影响。

#### Results
1. UMAP图并未表现出复杂的流型，点阵较为圆润，和教程上的UMAP图存在差异。这意味着UMAP结果不佳？
2. 文件夹`umap_clustering_figures_123`中的图片显示，批次效应不明显，可能是UMAP算法效果不好隐藏了批次效应。
3. `sc.pp.regress_out`在scanpy官方文档上的描述为 _Regress out (mostly) unwanted sources of variation. Uses simple linear regression. This is inspired by Seurat’s regressOut function in R. Note that this function tends to overcorrect in certain circumstances as described in issue526._ 说明还有可能是**批次效应被`sc.pp.regress_out(pp_adata, keys=["total_counts"])`解决了。虽然没有把keys指定为`batches`，但是`total_counts`确实也是三个批次存在细微差别的部分之一，可能起到了消除批次效应的作用。**
4. 老师注释的细胞类型准确，分割较为明显。有部分细胞类型占比极小。
5. UMAP图并没有展示出明显的时序变化。

### 脚本`insight_sc_data.py`
#### 内容
因为批次和n_genes_by_counts高度相关，编写了一个脚本观察批次、细胞发育时段、n_genes_by_n_counts三者的相关关系。
以及其它的绘图代码，脚本产生的所有图片均保存在脚本中的`FIG_FOLDER_PATH`文件夹中.

#### Results
1. 这三者高度相关。
2. 可能批次会影响n_genes_by_counts，有些批次测量很敏感，能测出很多种基因。
3. 可能细胞存在的发育时段也会影响n_genes_by_counts，细胞发育到某个时期表达的基因种类居然变多。
4. 批次三的所有细胞均来自E70时期。
5. 批次一和批次二的细胞在五个发育时间点之间分布得较为均匀

### 脚本`annotate_cell.py`
#### 内容
细胞类型注释脚本。

# 关于`sc.pp.scale`的理解和debug

我找到问题了！现在是2025年3月5日的晚上，找到scale函数的问题是因为`sc.tl.neighbors`函数构建的knn邻接图中存在大量的边权重恰好为1.0的边，而在图中可以明显观察到均值在0.22左右的钟形曲线，说明这些边均为异常值。
全程摸排数据预处理过程中的操作，最终发现`sc.pp.scale`函数存在异常。

首先让我们阅读一下sc.pp.scale的函数签名：
`scanpy.pp.scale(data, *, zero_center=True, max_value=None, copy=False, layer=None, obsm=None, mask_obs=None)`
这个函数对**基因（列）**进行标准化，使得**每个基因在各单细胞中的表达量的均值为0，方差为1**.
其中有一个参数叫max_value，这个参数的描述是这样子的：
```
max_value
float | None (default: None)
Clip (truncate) to this value after scaling. If None, do not clip.
```

并且需要注意的是，这个函数没有参数叫min_value，所以max_value很容易被理解为绝对值，即在经过scale后，数据的均值为0，方差为1后，那些特别大的值，大于max_value的值，会被截断到max_value。因为没有一个值叫min_value，因此很容易把max_value理解为绝对值，即那些特别小的值，小于-max_value的值，会被截断到-max_value。

然而实际上是错误的，对于小于-max_value的值，并不存在这样的截断机制。为什么我会得到这样的结论？
我在运行`sc.pp.scale`前后均执行了两次`print`操作：
```
    print(pd.Series(pp_adata.X.sum(axis=1).flatten()).describe())
    print(pd.Series(pp_adata.X.sum(axis=0).flatten()).describe())
```
这两次操作的意思是，分别对Anndata按行求和和按列求和，并且对求和后得到的一维数组调用`pandas.describe`方法计算统计量。
很明显，因为函数对**基因（列）**进行标准化，按行求和后计算得到的一维数组的统计量`print(pd.Series(pp_adata.X.sum(axis=0).flatten()).describe())`的输出很重要。
我的数据集中，当预处理到scale这一步时，基因个数有14750个，细胞个数有13605个。
当`sc.pp.scale(pp_adata, max_value=10)`时，缩放前的输出为：
```
count    13605.000000
mean        -0.002499
std         85.503703
min       -503.267315
25%        -55.706922
50%          0.669205
75%         56.809686
max        368.334427
dtype: float64
count    14749.000000
mean        -0.002305
std          0.093003
min         -3.318467
25%         -0.004804
50%          0.000395
75%          0.006740
max          2.227235
dtype: float64
```
缩放后的输出为：
```
count    13605.000000
mean       -68.099010
std        352.627792
min      -1437.517181
25%       -300.140177
50%        -82.132927
75%        146.897955
max       1782.512804
dtype: float64
count    1.474900e+04
mean    -6.281694e+01
std      1.073742e+02
min     -3.300094e+02
25%     -7.439085e+01
50%     -4.547474e-13
75%      2.486900e-14
max      2.241953e+01
dtype: float64
```

我们可以惊讶的发现，在缩放后，按行求和得到的一维数组的统计量非常奇怪。**均值为-60，最小值为-300。上四分位数非常接近0，而最大值居然为22？**
数组内的所有元素本应该非常趋近于0，而目前的情况并非如此，说明scale函数的结果存在很大问题。

当`sc.pp.scale(pp_adata, max_value=100)`时，缩放前的输出和上文展示的缩放前的输出一模一样。
而缩放后的输出：
```
count    1.360500e+04
mean     5.348003e-15
std      3.630238e+02
min     -1.341649e+03
25%     -2.431933e+02
50%     -1.857274e+01
75%      2.187736e+02
max      1.951758e+03
dtype: float64
count    1.474900e+04
mean     2.243962e-15
std      1.925448e-13
min     -8.668621e-13
25%     -8.881784e-14
50%      0.000000e+00
75%      9.592327e-14
max      8.135714e-13
dtype: float64
```
均值和方差和最大值和最小值和上下四分位数都非常接近0，说明该数组的所有元素都非常趋近于0。说明scale可能在正常运行了。

当`sc.pp.scale(pp_adata, max_value=None)`时，缩放前的输出和上文展示的缩放前的输出一模一样。
而缩放后的输出也和上文展示的一模一样，说明`max_value=100`足够大了。

由此我判断，对于小于`-max_value`的值，并不存在这样的截断机制，max_value只负责标准化后特别大的数。

经过上面的控制变量试验后，我发现对于这个数据集，`max_value`参数的设置对`sc.pp.scale`函数的影响极其巨大。

当把`max_value`参数拿掉之后，函数`sc.tl.neighbors`的异常未消失，只是图存在些许扰动，看来邻接图的异常不是由于`sc.pp.scale`导致的。

可以看看在计算leiden时抛出的exception总数，和异常边的总数是否一样。