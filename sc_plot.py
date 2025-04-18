import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
from typing import Optional, Callable, Tuple
from anndata import AnnData
from typing import Union
from os import PathLike

from matplotlib.pyplot import title


# ALL the function below should be used in the 'if __name__ == "__main__"' of sc_plot.py


def save_figure(filename: Union[str, PathLike]):
    """
    Saving figures to FIG_FOLDER_PATH/filename. This function is used in ALL remaining polt functions in this script.

    :param filename: given figure name, should include file extension.
    :return: None
    """
    if not os.path.exists(FIG_FOLDER_PATH):
        os.makedirs(FIG_FOLDER_PATH)  # 如果文件夹不存在，则创建
    filepath = os.path.join(FIG_FOLDER_PATH, filename)
    plt.savefig(filepath, bbox_inches="tight", dpi=300)  # 保存图片
    plt.close()  # 关闭当前图形，避免内存泄漏


def plot_scatter(adata: AnnData, x, y, title=None,
                 filename: Union[str, PathLike] = "scatter_n_genes_by_total_counts.png"):
    """
    Generate a scatter plot visualizing two observation annotations using scanpy.pl.scatter.

    :param adata: single cell sequencing data, anndata.AnnData object.
    :param x: Observation metadata column name for x-axis coordinates. Must exist in adata.obs.
    :param y: Observation metadata column name for y-axis coordinates. Must exist in adata.obs.
    :param title: Descriptive title displayed at the top of the plot. If None, no title will be shown.
    :param filename: Saving figure filename, file extension should include.
    :return: None
    """
    # scanpy中的绘图函数的保存逻辑和matplotlib不一样，虽然依赖于matplotlib
    sc.pl.scatter(adata, x=x, y=y, title=title, show=False, color="batch")
    save_figure(filename)


def plot_cell_counts_scatterplot(adata: AnnData, y_key="batch", x_key="cell_time",
                                 filename: Union[str, PathLike] = "batch_counts.png",
                                 is_proportion=False):
    """
    Generate a scatter plot visualizing cell numbers/proportion in
    two dimension of observation annotation using seaborn.scatterplot.

    :param is_proportion: Set True to plot cell proportion, set False to plot cell numbers.
    :param adata: AnnData object
    :param x_key: str, Observation metadata column name for x-axis coordinates. Must exist in adata.obs.
    :param y_key: str, Observation metadata column name for y-axis coordinates. Must exist in adata.obs.
    :param filename: str, Saving figure filename, should include file extension.
    :return: None (saves plot to file)
    """
    if y_key not in adata.obs_keys() or x_key not in adata.obs_keys():
        raise ValueError(f"'{y_key}' 或 '{x_key}' 不是adata对observation的注释。请检查前序工作是否完成")

    # 计算每个批次中每个细胞发育阶段的单细胞数量
    counts = adata.obs.groupby([x_key, y_key]).size().unstack(fill_value=0)
    if x_key == "cell_time":
        counts = counts.reindex(["E6.5", "E6.75", "E7.0", "E7.25", "E7.5"])

    if is_proportion:
        counts = counts.div(counts.sum(axis=1), axis=0)

    # 将计数数据转换为长格式，方便绘图
    counts_long = counts.reset_index().melt(id_vars=x_key, var_name=y_key, value_name="count")

    # Filter zeros out
    counts_long = counts_long[counts_long["count"] > 0]

    if y_key == "cell_time":
        counts_long[y_key] = pd.Categorical(
            counts_long[y_key],
            categories=["E6.5", "E6.75", "E7.0", "E7.25", "E7.5"],
            ordered=True
        )

    if x_key == "cell_type":
        counts_long[x_key] = pd.Categorical(
            counts_long[x_key],
            categories=["Epiblast", "Primitive Streak", "Nascent mesoderm", "Mixed mesoderm", "Mesenchyme",
                        "Haematoendothelial progenitors", "Blood progenitors 1", "Blood progenitors 2",
                        "Anterior Primitive Streak", "Def. endoderm", "Gut", "Visceral endoderm", "ExE endoderm"],
            ordered=True
        )

    # 使用 seaborn 绘制点图
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=counts_long,
        x=x_key,
        y=y_key,
        size="count",
        sizes=(100, 1000),  # 调整点的大小范围
        hue="count",  # 根据数量设置颜色
        palette="viridis",  # 颜色方案
        legend="full"
    )

    # 添加标题和标签
    plt.title(f"Batch Count by {x_key}")
    plt.xticks(rotation=45, ha='right', va='top', rotation_mode='anchor')
    plt.xlabel(x_key)
    plt.ylabel(y_key)
    title = "Count"
    if is_proportion:
        title = "Proportion"
    plt.legend(title=title, bbox_to_anchor=(1.05, 1), loc='upper left')

    # 保存图片
    save_figure(filename)


def plot_violin_n_genes_by_batch(adata: AnnData, filename: Union[str, PathLike] = "violin_by_batch.png"):
    """
    Plot violin polt of sequencing depth in each batch, showing the distribution of sequencing depth.

    :param adata: single cell sequencing data, anndata.AnnData object.
    :param filename: str, Saving figure filename, should include file extension.
    :return: None (saves plot to file)
    """
    if "batch" not in adata.obs_keys() or "n_genes_by_counts" not in adata.obs_keys():
        raise ValueError(f"batch 或 n_genes_by_counts 不是adata对observation的注释。请检查前序工作是否完成")

    # 创建子图，横向排列三张小提琴图
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 10), sharey=True)
    batches = adata.obs["batch"].unique()

    # 遍历每个批次，绘制小提琴图并输出前 5% 最大值
    for i, batch in enumerate(batches):
        # 筛选当前批次的细胞数据
        subset = adata.obs[adata.obs["batch"] == batch]

        # 计算前 5% 最大值（95% 分位数）
        top_5_percent_value = np.percentile(subset["n_genes_by_counts"], q=95)
        top_10_percent_value = np.percentile(subset["n_genes_by_counts"], q=90)
        print(f"Batch {batch}: 前 5% 最大值 = {top_5_percent_value}")
        print(f"Batch {batch}: 前 10% 最大值 = {top_10_percent_value}")

    # 遍历每个批次，绘制小提琴图
    for i, batch in enumerate(batches):
        # 筛选当前发育时段的细胞数据
        subset = adata.obs[adata.obs["batch"] == batch]

        # 绘制小提琴图
        sns.violinplot(
            y=subset["n_genes_by_counts"],
            ax=axes[i],
            color="skyblue",  # 设置颜色
            inner="quartile"  # 显示四分位数
        )

        # 设置标题和标签
        axes[i].set_title(f"Batch: {batch}")
        axes[i].set_xlabel(f"Batch: {batch}")
        axes[i].set_ylabel("n_genes_by_counts" if i == 0 else "")  # 只在第一个子图显示 y 轴标签

    # 调整布局
    plt.tight_layout()

    # 保存图片
    save_figure(filename)


def plot_grouped_violin(
        adata: AnnData,
        group_key: str,
        value_key: Optional[str] = None,
        value_calculator: Optional[Callable[[AnnData], pd.Series]] = None,
        title_prefix: str = "",
        filename: Union[str, PathLike] = "violin_plot.png",
        figsize: Tuple[int, int] = (15, 6),
        color: str = "skyblue",
        clip_quantile: Optional[float] = None,
        **kwargs
) -> None:
    """
    Generate grouped violin plots for observation metrics or calculated values.

    Can handle both precomputed observation metrics (value_key) and dynamically
    calculated values (value_calculator) across specified groups.

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix with observations and variables
    group_key : str
        Column name in adata.obs for grouping cells, such as cell type, cell time.
    value_key : str, optional
        Precomputed numeric column in adata.obs to visualize, such as n_genes_by_counts(sequencing depth)
    value_calculator : Callable, optional
        Function that takes subset AnnData and returns pd.Series of values
    title_prefix : str
        Text prefix for subplot titles (e.g., "Batch", "Cell Time")
    filename : str
        Output filename with extension (e.g., .png/.pdf)
    figsize : Tuple[int, int]
        Base dimensions for figure layout
    color : str
        Matplotlib color name for violins
    clip_quantile : float
        Percentile (0-100) for upper value truncation
    **kwargs
        Additional arguments passed to sns.violinplot

    Returns:
    --------
    None (saves plot to file)

    Examples:
    ---------
    # Plot precomputed observation metric
    plot_grouped_violin(adata, group_key='cell_time',
                       value_key='n_genes_by_counts',
                       title_prefix='Cell Time')

    # Plot calculated gene frequencies
    plot_grouped_violin(adata, group_key='batch',
                       value_calculator=lambda x: pd.Series(x.X.sum(axis=0).A1),
                       clip_quantile=95,
                       title_prefix='Batch')
    """
    # Validate input parameters
    if group_key not in adata.obs:
        raise ValueError(f"Group key '{group_key}' not found in adata.obs")

    if (value_key, value_calculator).count(None) != 1:
        raise ValueError("Must specify exactly one of: value_key or value_calculator")

    if value_key and value_key not in adata.obs:
        raise ValueError(f"Value key '{value_key}' not found in adata.obs")

    # Get unique groups and determine plot layout
    groups = adata.obs[group_key].unique()
    n_groups = len(groups)
    fig, axes = plt.subplots(
        nrows=1,
        ncols=n_groups,
        figsize=(figsize[0] * max(n_groups, 1), figsize[1]),
        sharey=True
    )
    axes = np.array(axes).flatten()  # Handle single-subplot case

    # Generate plots for each group
    for idx, (group, ax) in enumerate(zip(groups, axes)):
        # Subset data
        mask = adata.obs[group_key] == group

        subset = adata[mask]

        # Calculate values
        if value_calculator:
            values = value_calculator(subset)
        else:
            values = subset.obs[value_key]

        # Apply value clipping
        if clip_quantile:
            cap = np.percentile(values, clip_quantile)
            values = np.clip(values, None, cap)

        # Create violin plot
        sns.violinplot(
            y=values,
            ax=ax,
            color=color,
            inner="quartile",
            **kwargs
        )

        # Configure subplot
        ax.set_title(f"{title_prefix}: {group}", pad=12)
        ax.set_xlabel("")
        # only label the y-axis on the most left graph.
        if idx == 0:
            y_label = value_key if value_key else "Computed Values"
            ax.set_ylabel(y_label, labelpad=10)

    # Finalize layout and save
    plt.tight_layout()
    save_figure(filename)
    plt.close()


def plot_grouped_violin_gene_frequency(adata: AnnData, groupby='cell_time',
                                       filename: Union[str, PathLike] = "gene_frequency_distribution_by_time.png"):
    """
    绘制groupby中基因表达频数的小提琴图。

    参数:
    :param filename: 保存图片的文件名。
    :param groupby: 分类方法，对每一列绘制一张小提琴图。
    :param adata: AnnData 对象，包含单细胞测序数据。

    """
    import numpy as np
    # 检查批次信息是否存在
    if groupby not in adata.obs:
        raise ValueError(f"'{groupby}' 不在 adata.obs 中。请提供存在的对行注释名。")

    # 获取批次列表
    value_list = adata.obs[groupby].unique()

    # 计算每个基因在每个批次中的表达频数
    gene_freq_dict = {value: [] for value in value_list}

    # 创建子图，横向排列张小提琴图
    fig, axes = plt.subplots(nrows=1, ncols=len(value_list), figsize=(15, 10), sharey=True)

    for value in value_list:
        # 获取当前批次的细胞数据
        time_data = adata[adata.obs[groupby] == value]

        # 计算每个基因在当前批次中的表达频数
        gene_freq = np.sum(time_data.X, axis=0)
        gene_freq_dict[value] = gene_freq

    # 遍历每个batch，绘制小提琴图
    for i, value in enumerate(value_list):
        # 将基因表达频数转换为 pandas.Series
        gene_freq_series = pd.Series(gene_freq_dict[value])

        # 截断极端高值（取 99% 分位数作为上限）
        upper_limit = np.percentile(gene_freq_series, q=95)
        gene_freq_series = np.clip(gene_freq_series, a_min=None, a_max=upper_limit)

        # 绘制小提琴图
        sns.violinplot(
            y=gene_freq_series,
            ax=axes[i],
            color="skyblue",  # 设置颜色
            inner="quartile"  # 显示四分位数
        )

        # 设置标题和标签
        axes[i].set_title(f"{groupby}: {value}")
        axes[i].set_xlabel(f"{groupby}: {value}")
        axes[i].set_ylabel("Gene Frequency" if value == value_list[0] else "")  # 只在第一个子图显示 y 轴标签

    plt.tight_layout()
    save_figure(filename)


def plot_pairplot(adata: AnnData, filename: Union[str, PathLike] = "pairplot.png"):
    """
    plot pair plot using seaborn.pairplot.

    :param adata: single cell sequencing data, anndata.AnnData object.
    :param filename: str, Saving figure filename, should include file extension.
    :return: None (saves plot to file)
    """
    if "batch" not in adata.obs_keys() or "cell_time" not in adata.obs_keys() or "n_genes_by_counts" not in adata.obs_keys():
        raise ValueError("'batch', 'cell_time' 或 'n_genes_by_counts' 不在 adata.obs 中")

    # 绘制散点图矩阵
    sns.pairplot(
        adata.obs[["batch", "cell_time", "n_genes_by_counts"]],
        hue="batch",  # 按批次分组
        palette="Set2",  # 颜色方案
        diag_kind="kde"  # 对角线显示核密度估计图
    )

    # 保存图片
    save_figure(filename)


def plot_pair_violin(adata: AnnData, filename: Union[str, PathLike] = "grouped_violin.png"):
    """
    plot pair plot using seaborn.pairplot.

    :param adata: single cell sequencing data, anndata.AnnData object.
    :param filename: str, Saving figure filename, should include file extension.
    :return: None (saves plot to file)
    """

    # if "batch" not in adata.obs_keys() or "cell_time" not in adata.obs_keys() or "n_genes_by_counts" not in adata.obs_keys():
    #     raise ValueError("'batch', 'cell_time' 或 'n_genes_by_counts' 不在 adata.obs 中")
    import numpy as np

    adata.obs["gene_frequency"] = np.sum(adata.X, axis=1)

    # 创建分组小提琴图
    plt.figure(figsize=(12, 6))
    sns.violinplot(
        x="cell_time",  # 横轴为细胞发育时段
        # y="n_genes_by_counts",  # 纵轴为 n_genes_by_counts
        y="gene_frequency",  # 纵轴为 gene_frequency
        hue="batch",  # 按批次分组
        data=adata.obs,  # 数据来源
        palette="Set3",  # 颜色方案
        inner="quartile",  # 显示四分位数
        split=True  # 将不同批次的小提琴图分开显示
    )

    # 添加标题和标签
    plt.title("Gene Frequency by Cell Time and Batch")
    plt.xlabel("Cell Time")
    plt.ylabel("Gene Frequency")
    plt.legend(title="Batch", bbox_to_anchor=(1.05, 1), loc='upper left')

    # 保存图片
    save_figure(filename)


def plot_gene_regulation_heatmap(adata: AnnData, feature_key: str, feature_values: list,
                                 save_csv_path="gene_regulation_heatmap_data.csv",
                                 pval_threshold: float = 0.05, figsize=(10, 8), cmap='bwr',
                                 save_fig_path: Union[str, PathLike] = "gene_regulation_heatmap.png"):
    """
    Plot a heatmap of up-regulated and down-regulated genes based on any cell feature.

    Parameters:
    - adata: AnnData object containing single-cell RNA-seq data and cell annotations.
    - feature_key: str, the key in adata.obs that specifies the cell feature (e.g., developmental time, cell type).
    - feature_values: list of floats, the values of cell feature. The sequence will the sequence displayed on the heatmap.
    - pval_threshold: float, significance threshold for identifying up/down-regulated genes.
    - figsize: tuple, size of the heatmap figure.
    - cmap: str, color map for the heatmap (default is 'bwr' for blue-white-red).

    Returns:
    - None (displays the heatmap).
    """
    import matplotlib.pyplot as plt
    from scipy.stats import ranksums
    import concurrent.futures
    # Check feature_values are valid
    for feature in feature_values:
        if feature not in adata.obs[feature_key].values.unique():
            raise KeyError(f"Feature '{feature}' not found in adata.var['{feature_key}'].")
    if len(feature_values) == 0:
        raise ValueError("No feature values provided.")
    if len(feature_values) < len(adata.obs[feature_key].values.unique()):
        print(f"Warning: Provided {feature_values} are not complete."
              f"There are {len(adata.obs[feature_key].values.unique())} unique {feature_key} annotation values for the observation.")
    if len(feature_values) > len(adata.obs[feature_key].values.unique()):
        raise ValueError(f"Provided {feature_values} overflow."
                         f"There are {len(adata.obs[feature_key].values.unique())} unique {feature_key} annotation values for the observation.")

    # ------------------------------------

    # 以下是没有使用多线程的代码
    # Initialize heatmap matrices
    # heatmap_matrix[i, j]: Number of up-regulated genes when comparing feature_values[i] vs feature_values[j]
    # heatmap_matrix[j, i]: Number of down-regulated genes when comparing feature_values[i] vs feature_values[j]
    # n_features = len(feature_values)
    # heatmap_matrix = np.zeros((n_features, n_features))  # Matrix for up-regulated genes
    #
    # # Iterate over all pairs of feature values
    # for i in range(0, n_features):
    #     for j in range(i+1, n_features):
    #         # Get cells for the two feature values
    #         cells_type_i = adata[adata.obs[feature_key] == feature_values[i], :]
    #         cells_type_j = adata[adata.obs[feature_key] == feature_values[j], :]
    #
    #         for gene in adata.var_names:
    #             # Perform ranksums test to assess significance
    #             _, pval = ranksums(cells_type_i[:, gene].X, cells_type_j[:, gene].X)
    #             if pval < pval_threshold:
    #                 # Calculate mean difference
    #                 mean_diff = np.mean(cells_type_i[:, gene].X) - np.mean(cells_type_j[:, gene].X)
    #                 if mean_diff > 0:
    #                     heatmap_matrix[i, j] = heatmap_matrix[i, j] + 1
    #                 else:
    #                     heatmap_matrix[j, i] = heatmap_matrix[j, i] - 1

    # ------------------------------------

    n_features = len(feature_values)
    heatmap_matrix = np.zeros((n_features, n_features))  # Matrix for up-regulated genes

    # Define a function to process each (i, j) pair
    def process_pair(i, j):
        cells_type_i = adata[adata.obs[feature_key] == feature_values[i], :]
        cells_type_j = adata[adata.obs[feature_key] == feature_values[j], :]
        up_count = 0
        down_count = 0

        for gene in adata.var_names:
            # Perform ranksums test to assess significance
            _, pval = ranksums(cells_type_i[:, gene].X, cells_type_j[:, gene].X)
            if pval < pval_threshold:
                # Calculate mean difference
                mean_diff = np.mean(cells_type_i[:, gene].X) - np.mean(cells_type_j[:, gene].X)
                if mean_diff > 0:
                    up_count += 1
                else:
                    down_count += 1

        return i, j, up_count, down_count

    # Use ThreadPoolExecutor to parallelize the outer loops
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for i in range(n_features):
            for j in range(i + 1, n_features):
                futures.append(executor.submit(process_pair, i, j))

        # Collect results and update heatmap_matrix
        for future in concurrent.futures.as_completed(futures):
            i, j, up_count, down_count = future.result()
            heatmap_matrix[i, j] += up_count
            heatmap_matrix[j, i] -= down_count

    # ------------------------------------

    # Plot the heatmap
    plt.figure(figsize=figsize)
    color_bar_max = np.max(np.abs(heatmap_matrix))
    plt.imshow(heatmap_matrix, cmap=cmap, vmin=(-1) * color_bar_max, vmax=color_bar_max)
    plt.colorbar(label='Number of Up/Down-Regulated Genes')

    # Set axis labels
    plt.xticks(np.arange(n_features), feature_values, rotation=45)
    plt.yticks(np.arange(n_features), feature_values)
    plt.xlabel(f'Feature Value (j): {feature_key}')
    plt.ylabel(f'Feature Value (i): {feature_key}')
    plt.title(f'Up/Down-Regulated Genes Heatmap for {feature_key}')

    # Add grid lines and annotations
    for i in range(0, n_features):
        for j in range(0, n_features):
            if i == j:
                continue
            if heatmap_matrix[i, j] != 0:
                plt.text(j, i, f'{int(heatmap_matrix[i, j])}', ha='center', va='center', color='black')

    save_figure(save_fig_path)
    np.savetxt(fname=save_csv_path, X=heatmap_matrix, fmt='%d')
    print(f"Fig successfully saved to{save_fig_path}.")


def plot_david_result_bar(up_csv_path: Union[str, PathLike], down_csv_path: Union[str, PathLike], value_key="%", up_n_annotation: int = 10, down_n_annotation: int = 10,
                          filename="david_result_barplot.png", f_title="Primitive Streak vs Nascent mesoderm"):
    def read_df(file_path: Union[str, PathLike], n_head):
        _, ext = os.path.splitext(file_path)
        if ext == ".csv":
            df: pd.DataFrame = pd.read_csv(file_path, header=0, delimiter="\t")
        elif ext == ".xlsx":
            df: pd.DataFrame = pd.read_excel(file_path, header=0)
        else:
            df: pd.DataFrame = pd.read_table(file_path, header=0, delimiter="\t")
        df = df.head(n_head)
        df.sort_values(by="%", ascending=False, inplace=True, axis=0)
        return df


    up_gene_df = read_df(up_csv_path, n_head=up_n_annotation)
    down_gene_df = read_df(down_csv_path, n_head=down_n_annotation)

    n_terms_up = len(up_gene_df.index)
    n_terms_down = len(down_gene_df.index)
    longest = 0

    for term in pd.concat([up_gene_df, down_gene_df])["Term"]:
        if len(term) > longest:
            longest = len(term)

    # 创建画布和子图
    fig = plt.figure(figsize=(10, 12))

    # 将画布分为四个大小相同的块
    # 上半部分：y 轴标签在左上块，图主体在右上块
    # (left, bottom, width, height)
    ax1_label = fig.add_axes((0.1, 0.40 + n_terms_down * 0.03, longest * 0.01, n_terms_up * 0.03))  # 左上块：y 轴标签
    ax1_plot = fig.add_axes((0.1 + longest * 0.01, 0.40 + n_terms_down * 0.03, longest * 0.01, n_terms_up * 0.03))   # 右上块：图主体

    # 下半部分：y 轴标签在右下块，图主体在左下块
    ax2_label = fig.add_axes((0.1 + longest * 0.01, 0.15, longest * 0.01, n_terms_down * 0.03))  # 右下块：y 轴标签
    ax2_plot = fig.add_axes((0.1, 0.15, longest * 0.01, n_terms_down * 0.03))   # 左下块：图主体

    if value_key == "PValue":
        up_gene_df[value_key] = -np.log10(up_gene_df[value_key])
        down_gene_df[value_key] = -np.log10(down_gene_df[value_key])

    # 绘制上半部分的图
    ax1_plot.barh(up_gene_df['Term'], up_gene_df[value_key], color='green', label='Up-regulated')
    ax1_plot.set_ylabel('Term (Up-regulated)', fontsize=18)  # 设置轴标签
    ax1_plot.legend(loc='lower right')
    ax1_plot.tick_params(axis='y', labelsize=18)  # 调大轴刻度标签字号
    ax1_label.axis('off')  # 隐藏左上块的轴

    # 绘制下半部分的图
    ax2_plot.barh(down_gene_df['Term'], down_gene_df[value_key], color='red', label='Down-regulated')
    ax2_plot.set_ylabel('Term (Down-regulated)', fontsize=18)  # 设置轴标签
    ax2_plot.yaxis.set_label_position("right")  # 将轴标签移动到右边
    ax2_plot.yaxis.tick_right()  # 将刻度标签移动到右边
    ax2_plot.legend(loc='lower left')
    ax2_plot.tick_params(axis='y', labelsize=18)  # 调大轴刻度标签字号
    ax2_label.axis('off')  # 隐藏右下块的轴

    # 设置标题和 x 轴标签
    ax1_plot.set_title(f"{f_title} Up Regulated", fontsize=18)
    ax2_plot.set_title(f"{f_title} Down Regulated", fontsize=18)
    ax1_plot.set_xlabel("-Log10 p-value", fontsize=18)
    ax2_plot.set_xlabel("-Log10 p-value", fontsize=18)

    # 保存图像
    save_figure(filename)


if __name__ == "__main__":
    # NOTE: ALL the function in this script should be used in the 'if __name__ == "__main__"' of sc_plot.py
    #       ALL the images produced will be saved in FIG_FOLDER_PATH via the function save_figure.

    FIG_FOLDER_PATH = "./insight_figures"  # 图片保存文件夹路径
    adata = sc.read_h5ad("h5ads/final_627_combat_re.h5ad")
    # adata.obs.rename(
    #     columns={"stage": "cell_time", "cell": "cell_ID", "sequencing.batch": "batch", "celltype": "cell_type",
    #              "sample": "embryo_ID"}, inplace=True)
    # 调用绘图函数并保存图片
    # plot_scatter(adata, x="total_counts", y="n_genes_by_counts", title="Total Counts vs n_genes_by_counts", filename="scatter_total_counts_vs_n_genes.png")
    # plot_cell_time_proportion(adata, filename="cell_time_proportion.png")
    # plot_violin_n_genes_by_batch(adata, filename="violin_by_batch.png")
    # plot_violin_by_time(adata, filename="violin_by_time.png")
    # plot_pairplot(adata, filename="pairplot.png")
    # plot_grouped_violin(adata, filename="grouped_violin.png")
    # plot_gene_frequency_by_batch(adata, batch_key="batch", filename="gene_frequency_distribution_by_batch_3.png")
    # plot_cell_time_counts(adata)
    # plot_gene_frequency_by_time(adata)
    # plot_grouped_violin(adata)
    # plot_cell_time_counts(adata)
    # plot_cell_counts_scatterplot(adata, y_key="batch", x_key="cell_time", filename="627_cell_time_batch_counts.png")
    # cell_type_lst_a = ["Epiblast", "Primitive Streak", "Anterior Primitive Streak",
    #                    "Def. endoderm", "Gut", "Visceral endoderm", "ExE endoderm"]
    # cell_type_lst_b = ["Epiblast", "Primitive Streak", "Nascent mesoderm", "Mixed mesoderm", "Mesenchyme",
    #                    "Haematoendothelial progenitors", "Blood progenitors 1", "Blood progenitors 2"]
    # cell_type_lst_c = list(set(cell_type_lst_a + cell_type_lst_b))
    # # NOTE: cell_50327 cell_79163 are outlier cells, cell_86931 filtered.
    # # NOTE: cell_29635 cell_51875 cell_50514 cell_837, filtered.
    # sub_adata = adata[adata.obs['cell_type'].isin(cell_type_lst_c), :]
    # outlier_cells = ["cell_50327", "cell_79163", "cell_46746", "cell_29635", "cell_51875", "cell_50514", "cell_837",
    #                  "cell_86931"]
    # sub_adata = sub_adata[~sub_adata.obs_names.isin(outlier_cells), :]
    #
    # plot_cell_counts_scatterplot(sub_adata, y_key="batch", x_key="cell_type", filename="627_cell_type_batch_counts.png")
    # plot_cell_counts_scatterplot(sub_adata, y_key="cell_time", x_key="cell_type",
    #                              filename="627_cell_time_type_counts.png")
    # plot_gene_regulation_heatmap(adata=adata, feature_key="cell_type", save_fig_path="type_before_correction.png",
    #                              save_csv_path="type_before_correction.csv",
    #                              feature_values=["Primitive_Streak", "Mixed_mesoderm", "Nascent_mesoderm",
    #                                                         "Anterior_Primitive_Streak", "Somitic_mesoderm",
    #                                                         "Pharyngeal_mesoderm", "Paraxial_mesoderm",
    #                                                         "Intermediate_mesoderm"])
    # plot_gene_regulation_heatmap(adata=adata, feature_key="cell_time", save_fig_path="time_before correction.png",
    #                              save_csv_path="time_before_correction.csv",
    #                              feature_values=["E65", "E675", "E70", "E725", "E75"])
    # plot_gene_regulation_heatmap(adata=adata, feature_key="embryo_ID", save_fig_path="embryo_before_correction.png",
    #                              save_csv_path="embryo_before_correction.csv",
    #                              feature_values=["1", "2", "3", "4", "5", "6", "7", "10", "14", "15",
    #                                              "18", "19", "20", "23", "26", "27", "30", "31", "32"],
    #                              figsize=(16, 16))
    # plot_gene_regulation_heatmap(adata=adata, feature_key="batch", save_fig_path="batch_before_correction.png",
    #                              save_csv_path="batch_before_correction.csv",
    #                              feature_values=["1", "2", "3"])
    plot_david_result_bar(up_csv_path="david_result/PSN_up.xlsx", down_csv_path="david_result/PSN_down.xlsx",
                          value_key="PValue", up_n_annotation=100, down_n_annotation=100, filename="new2.png",
                          f_title="Primitive Streak vs Anterior Primitive Streak")
