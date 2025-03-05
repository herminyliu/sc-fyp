import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

from networkx.readwrite.json_graph.adjacency import adjacency_data

import preprocess  # 假设 preprocess 是你的自定义模块

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


def plot_cell_time_proportion(adata, filename="cell_time_proportion.png"):
    """绘制细胞发育时段比例图并保存"""
    if "cell_time" not in adata.obs.columns or "batch" not in adata.obs.columns:
        raise ValueError("'cell_time' 或 'batch' 不在 adata.obs 中")

    # 计算每个批次中每个细胞发育时段的比例
    counts = adata.obs.groupby(["batch", "cell_time"]).size().unstack(fill_value=0)
    proportions = counts.div(counts.sum(axis=1), axis=0)

    # 将比例数据转换为长格式，方便绘图
    proportions_long = proportions.reset_index().melt(id_vars="batch", var_name="cell_time", value_name="proportion")

    # 使用 seaborn 绘制点图
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=proportions_long,
        x="batch",
        y="cell_time",
        size="proportion",
        sizes=(100, 1000),  # 调整点的大小范围
        hue="proportion",  # 根据比例设置颜色
        palette="viridis",  # 颜色方案
        legend="full"
    )

    # 添加标题和标签
    plt.title("Cell Time Proportion by Batch")
    plt.xlabel("Batch")
    plt.ylabel("Cell Time")
    plt.legend(title="Proportion", bbox_to_anchor=(1.05, 1), loc='upper left')

    # 保存图片
    save_figure(filename)


def plot_cell_time_counts(adata, filename="cell_time_counts.png"):
    """绘制细胞发育时段比例图并保存"""
    if "cell_time" not in adata.obs.columns or "batch" not in adata.obs.columns:
        raise ValueError("'cell_time' 或 'batch' 不在 adata.obs 中")

    # 计算每个批次中每个细胞发育阶段的单细胞数量
    counts = adata.obs.groupby(["batch", "cell_time"]).size().unstack(fill_value=0)

    # 将计数数据转换为长格式，方便绘图
    counts_long = counts.reset_index().melt(id_vars="batch", var_name="cell_time", value_name="count")

    # 使用 seaborn 绘制点图
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=counts_long,
        x="batch",
        y="cell_time",
        size="count",
        sizes=(100, 1000),  # 调整点的大小范围
        hue="count",  # 根据数量设置颜色
        palette="viridis",  # 颜色方案
        legend="full"
    )

    # 添加标题和标签
    plt.title("Cell Time Count by Batch")
    plt.xlabel("Batch")
    plt.ylabel("Cell Time")
    plt.legend(title="Count", bbox_to_anchor=(1.05, 1), loc='upper left')

    # 保存图片
    save_figure(filename)


def plot_violin_n_genes_by_batch(adata, filename="violin_by_batch.png"):
    """按批次绘制小提琴图并保存"""
    # 创建子图，横向排列三张小提琴图
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10, 10), sharey=True)
    batches = adata.obs["batch"].unique()

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


if __name__ == "__main__":
    # 读取数据
    # origin_adata = sc.read_text(preprocess.FILE_PATH, delimiter='\t')
    # transposed_X = origin_adata.X.T
    # adata = sc.AnnData(X=transposed_X, var=origin_adata.obs, obs=origin_adata.var)

    # 注释细胞信息
    # adata = preprocess.annotate_cells(adata)
    adata = sc.read_h5ad("data_annotated_123.h5ad")

    # adata.write_h5ad(filename="data_annotated_123.h5ad")

    # 计算质量控制指标
    # sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

    # 调用绘图函数并保存图片
    # plot_scatter(adata, x="total_counts", y="n_genes_by_counts", title="Total Counts vs n_genes_by_counts", filename="scatter_total_counts_vs_n_genes.png")
    # plot_cell_time_proportion(adata, filename="cell_time_proportion.png")
    # plot_violin_by_batch(adata, filename="violin_by_batch.png")
    # plot_violin_by_time(adata, filename="violin_by_time.png")
    # plot_pairplot(adata, filename="pairplot.png")
    # plot_grouped_violin(adata, filename="grouped_violin.png")
    # plot_gene_frequency_by_batch(adata, batch_key="batch", filename="gene_frequency_distribution_by_batch_3.png")
    # plot_cell_time_counts(adata)
    # plot_gene_frequency_by_time(adata)
    plot_grouped_violin(adata)

