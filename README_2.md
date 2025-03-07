这个README文件是DeepSeek整理版，结构会更加清晰些，用语也会更加理性。

# 文件说明

由于细胞数据存在明显的批次效应，27 号的分析将三个批次分开进行了预处理和聚类。

---

## 原始数据文件

- **原始数据文件**：`Expr_Mat_PS+Mesoderm_Filter_ComBat_No.Rep_Cell.Name.Merge_2.txt`
- **测试数据文件**：`test_data.txt`  
  - 仅包含原始数据文件中前 4 个细胞（均来自第一批次），用于代码测试。

---

## 结果文件

### 预处理和聚类结果
1. **预处理后的 h5ad 文件**：`preprocessed_data_{batch_ID}.h5ad`  
   - 分别对应三个批次的预处理结果。
2. **预处理图表文件夹**：`preprocessed_figures_{batch_ID}`  
   - 包含预处理过程中生成的图表。
3. **聚类后的 h5ad 文件**：`clustered_data_{batch_ID}.h5ad`  
   - 分别对应三个批次的聚类结果。
4. **聚类图表文件夹**：`umap_clustering_figures_{batch_ID}`  
   - 包含聚类过程中生成的图表。

### 批次效应消除前的结果
1. **合并预处理后的 h5ad 文件**：`preprocessed_data_123.h5ad`  
   - 三个批次直接合并后的预处理结果。
2. **合并预处理图表文件夹**：`preprocessed_figures_123`  
   - 包含合并预处理过程中生成的图表。
3. **合并聚类后的 h5ad 文件**：`clustered_data_123.h5ad`  
   - 三个批次直接合并后的聚类结果。
4. **合并聚类图表文件夹**：`umap_clustering_figures_123`  
   - 包含合并聚类过程中生成的图表。

### 其他结果
1. **`insight_figures` 文件夹**：  
   - 包含 `insight_n_genes_by_counts.py` 脚本生成的结果图片。
2. **`{}_noBatchEffect` 后缀**：  
   - 表示在 `preprocess.py` 中消除批次效应后的结果。

---

## h5ad 文件对 `obs` 和 `var` 的注释说明

### 对 `obs` 的注释
- `batch`：细胞测序批次，值为 **1、2、3**，用于区分不同批次的实验数据。
- `cell_type`：注释的细胞类别，用于标识每个细胞的类型。
- `cell_time`：细胞发育时期，取值为 **E65、E675、D70、E725、E75**。
- `cell_ID`：细胞的唯一标识符。
- `n_genes`：每个细胞中检测到的基因数量。
- `n_genes_by_counts`：基于表达量（counts）检测到的基因数量。
- `total_counts`：每个细胞的总表达量（counts）。
- `leiden`：基于 Leiden 聚类算法得到的细胞聚类标签。

### 对 `var` 的注释
- `n_cells`：表达该基因的细胞数量。
- `n_cells_by_counts`：基于表达量（counts）检测到该基因的细胞数量。
- `mean_counts`：该基因在所有细胞中的平均表达量。
- `pct_dropout_by_counts`：该基因在细胞中的丢失率（dropout rate）。
- `total_counts`：该基因在所有细胞中的总表达量。
- `highly_variable`：是否被标记为高变基因（True/False）。
- `means`：标准化后的平均表达量。
- `dispersions`：基因的离散度（dispersion）。
- `dispersions_norm`：归一化离散度（normalized dispersion）。
- `mean`：标准化后的平均表达量（与 `means` 相同）。
- `std`：标准化后的标准差。

---

## 脚本文件

### `preprocess.py`
#### 编写参考
按照 [Scanpy 标准数据处理流程](https://scanpy.readthedocs.io/en/1.11.x/tutorials/basics/clustering-2017.html#finding-marker-genes) 编写。

#### 注意事项
1. **函数顺序**：Scanpy 的工作流存在依赖关系，请勿随意调整函数顺序。
2. **PCA 的必要性**：虽然聚类和拟时间分析使用 UMAP，但 UMAP 需要先运行 PCA。
3. **数据质量**：`Expr_Mat_PS+Mesoderm_Filter_ComBat_No.Rep_Cell.Name.Merge_2.txt` 已通过质量控制。
4. **批次效应**：在 `umap_clustering_figures_123` 中观察到明显的批次效应，三批细胞的 `n_genes_by_counts` 分布差异显著。
5. **批次效应处理**：同时使用 `sc.pp.regress_out` 和 `sc.pp.combat` 去除批次效应。
6. **提取批次**：提取第一批次细胞时，请使用 `adata_1 = adata[adata.obs.index.str.split('.').str[-2] == '1', :]`。

---

### `umap_clustering.py`
#### 编写参考
按照 [Scanpy 标准数据处理流程](https://scanpy.readthedocs.io/en/1.11.x/tutorials/basics/clustering-2017.html#finding-marker-genes) 编写。

#### 内容
包含对单细胞数据进行降维的代码，使用 UMAP 算法。但 UMAP 降维效果不佳，未学习到明显的流型。

---

### `show_batch_effect.py`
用于展示批次效应可能带来的影响。

#### 结果
1. **UMAP 图**：未表现出复杂流型，点阵较为圆润，与教程中的 UMAP 图存在差异。
2. **批次效应**：`umap_clustering_figures_123` 中的图片显示批次效应不明显，可能是 UMAP 效果不佳隐藏了批次效应。
3. **`sc.pp.regress_out`**：可能通过 `total_counts` 间接消除了批次效应。
4. **细胞类型**：老师注释的细胞类型准确，但部分类型占比极小。
5. **时序变化**：UMAP 图未展示出明显的时序变化。

---

### `insight_sc_data.py`
#### 内容
观察批次、细胞发育时段、`n_genes_by_counts` 三者的相关关系。生成的图片保存在 `FIG_FOLDER_PATH` 文件夹中。

#### 结果
1. 三者高度相关。
2. 批次可能影响 `n_genes_by_counts`，某些批次测量更敏感。
3. 细胞发育时段可能影响 `n_genes_by_counts`。
4. 批次三的所有细胞均来自 E70 时期。
5. 批次一和批次二的细胞在五个发育时间点分布均匀。

---

### `annotate_cell.py`
#### 内容
细胞类型注释脚本。

---

## 关于 `sc.pp.scale` 的理解和 Debug

### 问题背景
在分析单细胞数据时，发现 `sc.tl.neighbors` 函数构建的 KNN 邻接图中存在大量边权重恰好为 1.0 的边。同时，图中均值在 0.22 左右的钟形曲线表明这些边为异常值。经过排查，最终发现问题可能与 `sc.pp.scale` 函数有关。

---

### `sc.pp.scale` 函数签名
```
scanpy.pp.scale(data, *, zero_center=True, max_value=None, copy=False, layer=None, obsm=None, mask_obs=None)
```

#### 关键参数说明
- **`zero_center`**：是否将数据均值中心化（默认 `True`）。
- **`max_value`**：缩放后截断的最大值。**注意**：该参数仅截断大于 `max_value` 的值，不截断小于 `-max_value` 的值。
- **`copy`**：是否返回数据的副本（默认 `False`）。

---

### 问题分析
在运行 `sc.pp.scale` 后，发现 KNN 邻接图中存在异常边权重。为了排查问题，进行了以下实验：

#### 实验设计
1. **数据准备**：在 `sc.pp.scale` 前后分别对 `Anndata` 对象按行和按列求和，并计算统计量。
2. **参数设置**：分别测试 `max_value=10`、`max_value=100` 和 `max_value=None` 的情况。
3. **统计量计算**：使用 `pandas.describe` 方法计算均值、标准差、最小值、最大值等统计量。

---

### 实验结果

#### 1. `max_value=10`
- **缩放前**：
  ```plaintext
  count    13605.000000
  mean        -0.002499
  std         85.503703
  min       -503.267315
  25%        -55.706922
  50%          0.669205
  75%         56.809686
  max        368.334427
  dtype: float64
  ```
- **缩放后**：
  ```plaintext
  count    13605.000000
  mean       -68.099010
  std        352.627792
  min      -1437.517181
  25%       -300.140177
  50%        -82.132927
  75%        146.897955
  max       1782.512804
  dtype: float64
  ```
  **问题**：缩放后按行求和的统计量异常，均值为 -60，最小值为 -300，说明 `sc.pp.scale` 的结果存在问题。

#### 2. `max_value=100`
- **缩放前**：与 `max_value=10` 的缩放前结果一致。
- **缩放后**：
  ```plaintext
  count    1.360500e+04
  mean     5.348003e-15
  std      3.630238e+02
  min     -1.341649e+03
  25%     -2.431933e+02
  50%     -1.857274e+01
  75%      2.187736e+02
  max      1.951758e+03
  dtype: float64
  ```
  **结论**：缩放后统计量接近 0，说明 `sc.pp.scale` 正常运行。

#### 3. `max_value=None`
- **缩放前**：与 `max_value=10` 的缩放前结果一致。
- **缩放后**：与 `max_value=100` 的缩放后结果一致。
  **结论**：`max_value=100` 足够大，无需进一步调整。

---

### 关键发现
1. **`max_value` 的影响**：`max_value` 参数对 `sc.pp.scale` 函数的影响巨大。当 `max_value` 设置过小时，缩放结果会出现异常。
2. **截断机制**：`max_value` 仅截断大于 `max_value` 的值，不截断小于 `-max_value` 的值。因此，`max_value` 并非绝对值的截断阈值。
3. **KNN 邻接图异常**：移除 `max_value` 后，`sc.tl.neighbors` 的异常未消失，说明邻接图异常并非由 `sc.pp.scale` 导致。

---

### 进一步验证
为了确认 `sc.pp.scale` 是否导致 KNN 邻接图异常，进行了以下验证：
1. **控制变量实验**：在 `sc.pp.scale` 前后分别运行 `sc.tl.neighbors`，观察异常边权重的变化。
2. **异常边统计**：计算 `sc.tl.neighbors` 中异常边的总数，并与 `sc.pp.scale` 的异常统计量进行对比。
3. **Leiden 聚类异常**：观察在计算 Leiden 聚类时抛出的异常数量，是否与异常边总数一致。

---

### 结论
1. **`sc.pp.scale` 的正确使用**：在单细胞数据分析中，`sc.pp.scale` 的 `max_value` 参数需要根据数据特性合理设置，避免因截断不当导致数据异常。
2. **KNN 邻接图异常来源**：KNN 邻接图中的异常边权重并非由 `sc.pp.scale` 直接导致，可能需要进一步排查其他预处理步骤或数据本身的问题。
3. **实验的重要性**：通过控制变量实验和统计量分析，可以有效定位和解决单细胞数据分析中的问题。

---

### 建议
1. **合理设置 `max_value`**：根据数据分布特性，选择合适的 `max_value` 值，避免缩放后数据异常。
2. **全面排查问题**：在发现异常时，应逐步排查预处理步骤，结合统计量和可视化结果定位问题。
3. **记录实验过程**：详细记录每次实验的参数设置和结果，便于后续分析和复现。

通过以上分析和实验，深入理解了 `sc.pp.scale` 的行为，并成功定位了 KNN 邻接图异常的原因。

### 源代码探究

为了了解`sc.pp.scale`的行为，最好的办法就是查看源代码。

#### **`sc.pp.scale` 函数中 `max_value` 参数的工作原理**

`sc.pp.scale` 是 `scanpy` 中用于数据标准化的函数，其核心功能是将数据缩放到单位方差和零均值（如果 `zero_center=True`）。`max_value` 参数用于在标准化后对数据进行截断，以防止极端值对后续分析的影响。以下是 `max_value` 参数的工作原理及其在函数调用链中的传递过程。

---

#### **1. `max_value` 参数的传递路径**

`max_value` 参数从 `sc.pp.scale` 函数开始，经过 `scale_array` 函数，最终传递到 `clip_array` 函数中。以下是具体的传递路径：

1. **`sc.pp.scale` 函数**：
   - `max_value` 参数作为 `sc.pp.scale` 的输入参数之一。
   - 在 `sc.pp.scale` 函数中，`max_value` 被传递给 `scale_array` 函数。

   ```
   return scale_array(
       data, zero_center=zero_center, max_value=max_value, copy=copy, mask_obs=mask_obs
   )
   ```

2. **`scale_array` 函数**：
   - `scale_array` 函数接收 `max_value` 参数，并在标准化完成后调用 `clip_array` 函数进行截断。
   - 如果 `max_value` 不为 `None`，`scale_array` 会将 `max_value` 和 `zero_center` 参数传递给 `clip_array`。

   ```
   if max_value is not None:
       X = clip_array(X, max_value=max_value, zero_center=zero_center)
   ```

3. **`clip_array` 函数**：
   - `clip_array` 函数接收 `max_value` 和 `zero_center` 参数，并根据这些参数对数据进行截断。
   - 如果 `zero_center=True`，截断范围是 `[-max_value, max_value]`；如果 `zero_center=False`，截断范围是 `[-inf, max_value]`。

   ```
   a_min, a_max = -max_value, max_value
   if X[r, c] > a_max:
       X[r, c] = a_max
   elif X[r, c] < a_min and zero_center:
       X[r, c] = a_min
   ```

---

#### **2. `max_value` 和 `zero_center` 的关系**

`max_value` 和 `zero_center` 参数共同决定了截断的范围和方式：

1. **`zero_center=True`**：
   - 数据会被均值中心化，因此截断范围是双向的，即 `[-max_value, max_value]`。
   - 例如：
     ```
     X = np.array([-3, -2, 0, 2, 3])
     X = clip_array(X, max_value=2, zero_center=True)
     # 结果： [-2, -2, 0, 2, 2]
     ```

2. **`zero_center=False`**：
   - 数据不会被均值中心化，因此截断范围是单向的，即 `[-inf, max_value]`。
   - 例如：
     ```
     X = np.array([-3, -2, 0, 2, 3])
     X = clip_array(X, max_value=2, zero_center=False)
     # 结果： [-3, -2, 0, 2, 2]
     ```

---

#### **3. 总结**

- **`max_value` 的作用**：在标准化后对数据进行截断，以防止极端值的影响。
- **传递路径**：`max_value` 从 `sc.pp.scale` 传递到 `scale_array`，最终传递到 `clip_array`。
- **与 `zero_center` 的关系**：
  - 如果 `zero_center=True`，截断范围是 `[-max_value, max_value]`。
  - 如果 `zero_center=False`，截断范围是 `[-inf, max_value]`。

通过这种设计，`sc.pp.scale` 函数可以灵活地处理不同场景下的数据标准化和截断需求。

