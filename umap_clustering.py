import scanpy as sc
from typing import Literal

FILE_PATH = "preprocessed_data_123_noBatchEffect.h5ad"
SAVING_FIG_FOLDER = './umap_clustering_figures_123_noBatchEffect'
SAVING_FILE_PATH = './clustered_data_123_noBatchEffect.h5ad'  # 应该为h5ad文件
N_MARKER_GENE = 25
IS_AUTOSAVE = False
IS_AUTOSHOW = True


def cluster_umap(u_adata, leiden_resolution):
    sc.pp.pca(u_adata, svd_solver=None, zero_center=True, n_comps=50)
    # sc.pl.pca_variance_ratio(u_adata, log=True, n_pcs=50)
    sc.pp.neighbors(u_adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(
        u_adata,
        resolution=leiden_resolution,
        random_state=42,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    sc.tl.paga(u_adata, groups="leiden")
    sc.pl.paga(u_adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
    # Plot PAGA first, so that `adata.uns['paga']['pos']` exists. Need to run sc.pl.paga before sc.pl.umap.
    sc.tl.umap(u_adata, init_pos='paga')
    return u_adata


def debug_pp_neighbors(u_adata):
    '''
    对sc.pp.neighbors(u_adata, n_neighbors, n_pcs)生成的KNN邻接矩阵进行表征
    :param u_adata: 传入的Anndata
    :return: None
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    sc.pp.pca(u_adata, svd_solver="arpack", n_comps=50)
    sc.pp.neighbors(u_adata, n_neighbors=10, n_pcs=40, metric='cosine')
    # 计算每个节点的度
    degrees = np.array(u_adata.obsp['connectivities'].sum(axis=1)).flatten()

    # 绘制度分布图
    plt.figure(figsize=(8, 6))
    plt.hist(degrees, bins=50, color='blue', alpha=0.7)
    plt.title('kNN Graph Degree Distribution')
    plt.xlabel('Degree')
    plt.ylabel('Frequency')
    plt.show()

    # 提取 kNN 图的权重
    weights = u_adata.obsp['connectivities'].data

    # 绘制权重分布图
    plt.figure(figsize=(8, 6))
    plt.hist(weights, bins=50, color='green', alpha=0.7)
    plt.title('kNN Graph Weight Distribution')
    plt.xlabel('Weight')
    plt.ylabel('Frequency')
    plt.show()

    # 提取稀疏矩阵
    connectivities = u_adata.obsp['connectivities']

    # 找到权重为 1.0 的边的行和列索引
    rows, cols = connectivities.nonzero()  # 获取所有非零元素的行和列索引
    weights = connectivities.data  # 获取所有非零元素的权重
    edges_with_weight_1 = (weights == 1.0)  # 找到权重为 1.0 的边

    # 统计每个节点的权重为 1.0 的边数
    node_degrees_weight_1 = np.zeros(connectivities.shape[0])  # 初始化每个节点的权重为 1.0 的边数
    for row, col in zip(rows[edges_with_weight_1], cols[edges_with_weight_1]):
        node_degrees_weight_1[row] += 1  # 统计每个节点的权重为 1.0 的边数

    # 绘制每个节点的权重为 1.0 的边数
    plt.figure(figsize=(8, 6))
    plt.scatter(range(len(node_degrees_weight_1)), node_degrees_weight_1, alpha=0.5)
    plt.title('Number of Weight=1.0 Edges per Node')
    plt.xlabel('Node Index')
    plt.ylabel('Number of Weight=1.0 Edges')
    plt.show()


def finding_marker_gene(m_adata, groupby, saving_fig_folder, n_marker_gene,
                        method: Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var"], pl_groups=None):
    sc.tl.rank_genes_groups(m_adata, groupby=groupby, groups="all", method=method)
    sc.pl.rank_genes_groups(m_adata, n_genes=n_marker_gene, sharey=False)
    sc.pl.rank_genes_groups_dotplot(
        m_adata, groupby=groupby, groups=pl_groups,
        standard_scale="var", n_genes=3,
    )
    from pandas import DataFrame
    if method == "t-test":
        DataFrame(m_adata.uns['rank_genes_groups']['logfoldchanges']).to_csv(f"./csv/rank_gene_{groupby}_{method}_lfg.csv")
    DataFrame(m_adata.uns['rank_genes_groups']['pvals']).to_csv(f"./csv/rank_gene_{groupby}_{method}_pvals.csv")
    DataFrame(m_adata.uns['rank_genes_groups']['names']).to_csv(f"./csv/rank_gene_{groupby}_{method}_names.csv")
    # pts会不时地报错，很讨嫌
    # DataFrame(m_adata.uns['rank_genes_groups']['pts']).to_csv(f"./csv/rank_gene_{groupby}_{method}_pts.csv")
    return m_adata


def plot_marker_gene(pm_adata, marker_genes_lst):
    # var_names should be a valid subset of adata.var_names.
    for gene in marker_genes_lst:
        if gene not in pm_adata.var_names:
            raise ValueError("Gene {} do not exist in the file.".format(gene))
    sc.pl.dotplot(adata=pm_adata, var_names=marker_genes_lst, groupby="leiden", save=f"leiden_dendrogram.pdf", dendrogram=True)
    sc.pl.dotplot(adata=pm_adata, var_names=marker_genes_lst, groupby="cell_time", save=f"cell_time_dendrogram.pdf", dendrogram=True)
    sc.pl.dotplot(adata=pm_adata, var_names=marker_genes_lst, groupby="cell_type", save=f"cell_type_dendrogram.pdf", dendrogram=True)


def plot_umap(pp_adata):
    import matplotlib.pyplot as plt
    sc.pl.umap(pp_adata, color=["cell_type", "cell_time", "batch", "leiden"], wspace=0.7, hspace=0.25, ncols=2)

    plt.tight_layout()


def plot_stack_violin(pp_adata, marker_genes_lst):
    sc.pl.stacked_violin(adata=pp_adata, var_names=marker_genes_lst, groupby="cell_type", dendrogram=False)


if __name__ == "__main__":
    # 设置展示运行中会出现的信息
    sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
    # 打印运行环境
    sc.logging.print_header()
    # 设置图像参数，保存信息
    sc.settings.set_figure_params(dpi=150, dpi_save=300, facecolor="white")
    sc.settings.figdir = SAVING_FIG_FOLDER
    sc.settings.autosave = IS_AUTOSAVE
    sc.settings.autoshow = IS_AUTOSHOW

    adata = sc.read_h5ad(FILE_PATH)
    debug_pp_neighbors(adata)
    # adata = finding_marker_gene(adata, "logreg")
    # adata = finding_marker_gene(adata, "t-test")
    # adata = finding_marker_gene(adata, "wilcoxon")
    # adata.write_h5ad(SAVING_FILE_PATH)
    marker_genes_lst = [
    "Gm16120",
    "Gm12446",
    "Vmn1r18",
    "Abca1",
    "Gm42705",
    "Ccdc171",
    "Gm38037",
    "Camk1g",
    "D2Bwg1423e",
    "Vps29",
    "Bicdl1",
    "Camk1g",
    "Gm13442"
]
    # plot_marker_gene(adata, marker_genes_lst)
    # plot_prior_info(adata)
    # plot_stack_violin(adata, marker_genes_lst)

