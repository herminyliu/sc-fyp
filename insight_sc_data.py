import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

# ALL the function below should be used in the 'if __name__ == "__main__"' of insight_sc_data.py

# 定义全局变量，指定图片保存路径
FIG_FOLDER_PATH = "./insight_figures"  # 图片保存文件夹路径


def save_figure(filename):
    """保存当前图片到指定文件夹"""
    if not os.path.exists(FIG_FOLDER_PATH):
        os.makedirs(FIG_FOLDER_PATH)  # 如果文件夹不存在，则创建
    filepath = os.path.join(FIG_FOLDER_PATH, filename)
    plt.savefig(filepath, bbox_inches="tight", dpi=300)  # 保存图片
    plt.close()  # 关闭当前图形，避免内存泄漏


def plot_scatter(adata, x, y, title=None, filename="scatter_n_genes_by_total_counts.png"):
    """绘制散点图并保存"""
    # scanpy中的绘图函数的保存逻辑和matplotlib不一样，虽然依赖于matplotlib
    sc.pl.scatter(adata, x=x, y=y, title=title, show=False, color="batch")
    save_figure(filename)


def plot_cell_time_counts(adata, y_key="batch", x_key="cell_time", filename="batch_counts.png", is_proprotion=False):
    """绘制细胞发育时段比例图并保存

    参数:
    adata: AnnData 对象，包含单细胞数据
    batch_key: str, 指定表示批次的列名，默认为 'batch'
    cell_time_key: str, 指定表示细胞发育时段的列名，默认为 'cell_time'
    filename: str, 保存图片的文件名，默认为 'batch_counts.png'
    """
    if y_key not in adata.obs.columns or x_key not in adata.obs.columns:
        raise ValueError(f"'{y_key}' 或 '{x_key}' 不在 adata.obs 中")

    # 计算每个批次中每个细胞发育阶段的单细胞数量
    counts = adata.obs.groupby([x_key, y_key]).size().unstack(fill_value=0)

    if is_proprotion:
        counts = counts.div(counts.sum(axis=1), axis=0)

    # 将计数数据转换为长格式，方便绘图
    counts_long = counts.reset_index().melt(id_vars=x_key, var_name=y_key, value_name="count")

    # 过滤掉值为0的数据点
    counts_long = counts_long[counts_long["count"] > 0]

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
    plt.xlabel(x_key)
    plt.ylabel(y_key)
    title = "Count"
    if is_proprotion:
        title = "Proportion"
    plt.legend(title=title, bbox_to_anchor=(1.05, 1), loc='upper left')

    # 保存图片
    save_figure(filename)


def plot_violin_n_genes_by_batch(adata, filename="violin_by_batch.png"):
    """按批次绘制小提琴图并保存"""
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
    # save_figure(filename)


def plot_violin_n_genes_by_time(adata, filename="violin_by_time.png"):
    """按发育时间绘制小提琴图并保存"""
    if "cell_time" not in adata.obs.columns or "n_genes_by_counts" not in adata.obs.columns:
        raise ValueError("'cell_time' 或 'n_genes_by_counts' 不在 adata.obs 中")

    # 获取细胞的发育时段（假设有五个时间点）
    cell_times = adata.obs["cell_time"].unique()
    if len(cell_times) != 5:
        raise ValueError("'cell_time' 必须有五个时间点")

    # 创建子图，横向排列五张小提琴图
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(15, 10), sharey=True)
    cell_times = adata.obs["cell_time"].unique()

    # 遍历每个发育时段，绘制小提琴图
    for i, cell_time in enumerate(cell_times):
        # 筛选当前发育时段的细胞数据
        subset = adata.obs[adata.obs["cell_time"] == cell_time]

        # 绘制小提琴图
        sns.violinplot(
            y=subset["n_genes_by_counts"],
            ax=axes[i],
            color="skyblue",  # 设置颜色
            inner="quartile"  # 显示四分位数
        )

        # 设置标题和标签
        axes[i].set_title(f"Cell_Time: {cell_time}")
        axes[i].set_xlabel(f"Cell_Time: {cell_time}")
        axes[i].set_ylabel("n_genes_by_counts" if i == 0 else "")  # 只在第一个子图显示 y 轴标签

    # 调整布局
    plt.tight_layout()

    # 保存图片
    save_figure(filename)


def plot_gene_frequency_by_batch(adata, batch_key='batch', filename="gene_frequency_distribution_by_batch.png"):
    """
    绘制每个批次中基因表达频数的小提琴图。

    参数:
    - adata: AnnData 对象，包含单细胞测序数据。
    - batch_key: str, adata.obs 中存储批次信息的列名，默认为 'batch'。
    - figsize: tuple, 图像的大小，默认为 (10, 6)。
    """
    import numpy as np
    # 检查批次信息是否存在
    if batch_key not in adata.obs:
        raise ValueError(f"'{batch_key}' 不在 adata.obs 中。请提供正确的批次信息列名。")

    # 获取批次列表
    batch_list = adata.obs[batch_key].unique()

    # 计算每个基因在每个批次中的表达频数
    gene_freq_dict = {batch: [] for batch in batch_list}

    # 创建子图，横向排列张小提琴图
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 10), sharey=True)

    for batch in batch_list:
        # 获取当前批次的细胞数据
        batch_data = adata[adata.obs[batch_key] == batch]

        # 计算每个基因在当前批次中的表达频数
        gene_freq = np.sum(batch_data.X, axis=0)
        gene_freq_dict[batch] = gene_freq


    # 遍历每个batch，绘制小提琴图
    for i, batch in enumerate(batch_list):
        # 将基因表达频数转换为 pandas.Series
        gene_freq_series = pd.Series(gene_freq_dict[batch])
        # print(gene_freq_dict[batch][:10])  # 打印前 10 个基因的表达频数
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
        axes[i].set_title(f"Batch: {batch}")
        axes[i].set_xlabel(f"Batch: {batch}")
        axes[i].set_ylabel("Gene Frequency" if batch == 0 else "")  # 只在第一个子图显示 y 轴标签

    plt.tight_layout()
    save_figure(filename)


def plot_gene_frequency_by_time(adata, time_key='cell_time', filename="gene_frequency_distribution_by_time.png"):
    """
    绘制每个批次中基因表达频数的小提琴图。

    参数:
    - adata: AnnData 对象，包含单细胞测序数据。
    - batch_key: str, adata.obs 中存储批次信息的列名，默认为 'batch'。
    - figsize: tuple, 图像的大小，默认为 (10, 6)。
    """
    import numpy as np
    # 检查批次信息是否存在
    if time_key not in adata.obs:
        raise ValueError(f"'{time_key}' 不在 adata.obs 中。请提供正确的批次信息列名。")

    # 获取批次列表
    time_list = adata.obs[time_key].unique()

    # 计算每个基因在每个批次中的表达频数
    gene_freq_dict = {time: [] for time in time_list}

    # 创建子图，横向排列张小提琴图
    fig, axes = plt.subplots(nrows=1, ncols=5, figsize=(15, 10), sharey=True)

    for time in time_list:
        # 获取当前批次的细胞数据
        time_data = adata[adata.obs[time_key] == time]

        # 计算每个基因在当前批次中的表达频数
        gene_freq = np.sum(time_data.X, axis=0)
        gene_freq_dict[time] = gene_freq


    # 遍历每个batch，绘制小提琴图
    for i, time in enumerate(time_list):
        # 将基因表达频数转换为 pandas.Series
        gene_freq_series = pd.Series(gene_freq_dict[time])
        # print(gene_freq_dict[batch][:10])  # 打印前 10 个基因的表达频数
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
        axes[i].set_title(f"Time: {time}")
        axes[i].set_xlabel(f"Time: {time}")
        axes[i].set_ylabel("Gene Frequency" if time == "E65" else "")  # 只在第一个子图显示 y 轴标签

    plt.tight_layout()
    save_figure(filename)


def plot_pairplot(adata, filename="pairplot.png"):
    """绘制散点图矩阵并保存"""
    if "batch" not in adata.obs.columns or "cell_time" not in adata.obs.columns or "n_genes_by_counts" not in adata.obs.columns:
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


def plot_grouped_violin(adata, filename="grouped_violin.png"):
    """绘制分组小提琴图并保存"""
    # if "batch" not in adata.obs.columns or "cell_time" not in adata.obs.columns or "n_genes_by_counts" not in adata.obs.columns:
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


def plot_gene_regulation_heatmap(adata, feature_key: str, feature_values: list, save_csv_path="gene_regulation_heatmap_data.csv",
                                 pval_threshold: float=0.05, figsize=(10, 8), cmap='bwr', save_fig_path="gene_regulation_heatmap.png"):
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

    plt.savefig(save_fig_path, dpi=300)
    np.savetxt(fname=save_csv_path, X=heatmap_matrix, fmt='%d')
    print(f"Fig successfully saved to{save_fig_path}.")


if __name__ == "__main__":
    # ALL the function in this script should be used in the 'if __name__ == "__main__"' of insight_sc_data.py
    # ALL the images produced will be saved in FIG_FOLDER_PATH via the function save_figure.
    adata = sc.read_h5ad("raw_627d.h5ad")
    adata.obs.rename(
        columns={"stage": "cell_time", "cell": "cell_ID", "sequencing.batch": "batch", "celltype": "cell_type",
                 "sample": "embryo_ID"}, inplace=True)
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
    plot_cell_time_counts(adata, filename="627_cell_type_batch_counts.png")
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
