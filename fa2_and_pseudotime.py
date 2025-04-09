import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt
from pandas import DataFrame
import igraph

FILE_PATH = "clustered_data_123_noBatchEffect.h5ad"
IS_AUTOSAVE = True
IS_AUTOSHOW = False


def paga(p_adata, color, groups):
    """
        Mapping out the coarse-grained connectivity structures of complex manifolds
        groupby is not in the params list of function sc.tl.paga
        illustration of param groups：
        Key for categorical in adata.obs. You can pass your predefined groups by choosing any categorical annotation of observations.

        groups and groupby is not in the params list of function sc.pl.paga, since groups is inherited from sc.tl.paga.
        illustration of param labels in sc.pl.paga:
        The node labels. If None, this defaults to the group labels stored in the categorical for which paga() has been computed.
    """
    sc.tl.paga(p_adata, groups=groups)
    sc.pl.paga(
        layout="fa", random_state=100,
        adata=p_adata, color=color, labels=None, save=f"paga_{groups}.png", fontsize=8, node_size_scale=1.0, node_size_power=0.5,
        title=f"PAGA Abstract Graph"
    )
    return p_adata


def paga_scatter(ps_adata: AnnData, color):
    """
        sc.pl.draw_graph绘制精细的paga散点图，图元素包含很多散点，每个点就是一个单细胞
        sc.pl.paga绘制的是粗略的图
        sc.tl.draw_graph不包含groups和groupby这两种参数。groups参数自动继承自前面的tl.paga
        sc.pl.draw_graph不包含groups和groupby这两种参数
    """

    sc.tl.draw_graph(ps_adata, init_pos="paga")
    (sc.get.obs_df(adata=ps_adata, obsm_keys=[('X_draw_graph_fa', 0), ('X_draw_graph_fa', 1)])
     .to_csv("temp4/graph_coordinate.csv"))


    sc.pl.draw_graph(
        adata=ps_adata, color=color, legend_loc='right margin', wspace=2.0,
        save=f"_force_directed_graph_{str(color[0])}.png",
        ncols=4
    )
    return ps_adata


def plot_fa_color_pseudotime(p_adata, color: list[str] = None, fa_titie: str = None, umap_title: str = None):
    """
    Plot 2D ForceAtlas manifold map and UMAP map coloring with p_adata.obs['dpt_order_indices']

    Parameters:
        adata (AnnData): Annotated data matrix.
        :param umap_title:
        :param fa_titie:
        :param p_adata:
        :param color:
    Returns:
        No returns

    """
    if color is None:
        color = ['normalized_dpt']
    if 'dpt_order_indices' not in p_adata.obs_keys():
        raise KeyError(f"'dpt_order_indices' not found. Please compute dpt first.")
    tl_drawgraph_layout = ['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa']
    tl_drawgraph_key = [f'X_draw_graph_{i}' for i in tl_drawgraph_layout]
    key_missing_flag = True
    for i in tl_drawgraph_key:
        if i in p_adata.obsm:
            key_missing_flag = False
    if key_missing_flag:
        raise KeyError(f"'X_draw_graph_['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa']' not found. Please compute tl.draw_graph first.")

    if 'normalized_dpt' in color:
        # Normalize 'dpt_order_indices' to a range of [0, 1] for color mapping
        dpt_order_indices = p_adata.obs['dpt_order_indices']
        normalized_dpt = (dpt_order_indices - dpt_order_indices.min()) / (dpt_order_indices.max() - dpt_order_indices.min())
        # Add normalized pseudotime to p_adata.obs for coloring
        p_adata.obs['normalized_dpt'] = normalized_dpt

    # Plot the 2D manifold map with custom coloring
    sc.pl.draw_graph(
        adata=p_adata,
        color=color,  # Use normalized pseudotime values for coloring
        color_map="viridis",        # Apply the colormap
        # title=fa_titie,
        show=False,
        save=f"{color[3:6]}.png",
        ncols=3
    )

    sc.pl.umap(
        p_adata,
        color=color,  # Use normalized pseudotime values for coloring
        color_map="viridis",
        # title=umap_title,
        show=False,
        save=f"{color[3:6]}.png",
        ncols=3)


def diff(dm_adata, color: list[str]):
    """
    Do diffusion map analysis and plot graph.

    :param dm_adata: Anndata object
    :param color: mapping colors to observation annotation to distinguish.
    :return:
    """
    sc.tl.diffmap(dm_adata, n_comps=15)
    sc.pl.diffmap(dm_adata, color=color,
                  ncols=4, save=f"diffmap_{str(color[0])}.png")
    return dm_adata


def pseudotime(dm_adata: AnnData, n_branchings):
    """
    Do pseudotime calculation using scanpy.tl.dpt. if n_branchings > 1, plot sc.pl.dpt_groups_pseudotime.

    Note: Cited from the scanpy docs:
        ------------------
        dpt() requires running neighbors(), first. dpt() also requires to run diffmap() first.
        As previously, dpt() came with a default parameter of n_dcs=10 but diffmap() has a default parameter of n_comps=15,
        you need to pass n_comps=10 in diffmap() in order to exactly reproduce previous dpt() results.
        ------------------
        while in this repo, neighbors() is run in paga and umap clustering while diffmap() is run in fa2 previously,
        this script assert dm_adata contain the results already.

    :param dm_adata: AnnData object
    :param n_branchings: Number of branchings to detect. See details in scanpy.tl.dpt.
    :return: AnnData object with newly added observation annotation about dpt result. Save a fig if n_branchings > 1
    """
    from numpy import flatnonzero, argsort, array

    # if calculated pp.neighbors before, should contain annotation "distances" "connectivities" both.
    # if (not any(name.startswith("distances") for name in dm_adata.obs_keys())
    #         or (not any(name.startswith("connectivities") for name in dm_adata.obs_keys()))):
    #     raise ValueError("adata have not computed sc.pp.neighbors() before.")
    #
    # if "X_diffmap" not in dm_adata.obsm_keys():
    #     raise ValueError("adata have not computed sc.tl.diffmap() before.")

    # Check whether root time exist.
    # Set all the cell in the earliest stage to be the root cells,which is only very small amount.
    root_time = 'Epiblast'
    if root_time not in dm_adata.obs['cell_type'].unique():
        raise ValueError(f"{root_time} not found in dm_adata.obs['cell_type']. Cannot set root cell.")

    # set root cell, must do before run dpt.
    dm_adata.uns['iroot'] = flatnonzero(dm_adata.obs['cell_type'] == root_time)[0]

    # NOTE:
    #   If n_branchings==0, no field adata.obs['dpt_groups'] adata.obs['dpt_order_indices'] will be written
    #   sc.pl.dpt_timeseries requires adata.obs['dpt_order_indices'] adata.uns['dpt_changepoint'] to plot
    #   sc.pl.dpt_groups_pseudotime requires:
    #   adata.obs['dpt_order_indices'] adata.obs['dpt_groups'] adata.uns['dpt_changepoint']to plot

    sc.tl.dpt(dm_adata, n_dcs=15, n_branchings=n_branchings, copy=False)

    print("===scl.tl.dpt has been computed====")

    if 'dpt_pseudotime' not in dm_adata.obs_keys():
        raise KeyError(f"'dpt_pseudotime' not found. DPT may not have run correctly.")

    if n_branchings >= 1 and (('dpt_groups' not in dm_adata.obs_keys()) or ('dpt_order_indices' not in dm_adata.obs_keys())):
        raise KeyError("dm_adata.obs['dpt_groups'] or dm_adata.obs['dpt_order_indices'] is None. something is wrong in sc.tl.dpt.")

    if n_branchings >= 1:
        sc.pl.dpt_groups_pseudotime(dm_adata, color_map="viridis", save=".png")

    if n_branchings == 0:
        # construct obs["dpt_order_indices"] by hand
        dm_adata.obs["dpt_order_indices"] = argsort(dm_adata.obs["dpt_pseudotime"])
        # construct uns["dpt_changepoints"] by hand, to show data contain no branch thus no change point.
        dm_adata.uns["dpt_changepoints"] = array([])

    return dm_adata


def plot_pseudotime(dm_adata: AnnData, gene_list: list[str] = None, cell_sample_step: int = 100, color_map="viridis"):
    """
    plot the expression of gene in gene_list over pesudotime using sc.pl.dpt_timeseries.

    :param dm_adata: Anndata object
    :param gene_list: Genes to plot on the heatmap.
    :param cell_sample_step: int, set to 1 to use all cells.
    :param color_map: color map used on heatmap.
    :return: None(fig saved)
    """
    from numpy import argsort
    if cell_sample_step < 1:
        print("Warning: cell_sample_step should be no less than 1, automatically set to 1.")
        cell_sample_step = 1
    missing_genes = [gene for gene in gene_list if gene not in dm_adata.var_names]
    if missing_genes:
            raise ValueError(f"The following genes are not found in the dataset: {missing_genes}")

    # 筛选基因
    adata_filtered = dm_adata[:, gene_list].copy()

    # 采样细胞
    cell_indices = adata_filtered.obs["dpt_order_indices"].values[::cell_sample_step]
    adata_sampled = adata_filtered[cell_indices].copy()
    # Rearrange obs["dpt_order_indices"]
    adata_sampled.obs["dpt_order_indices"] = argsort(adata_sampled.obs["dpt_pseudotime"])
    sc.pl.dpt_timeseries(adata_sampled, color_map=color_map, save=".png", as_heatmap=True)


def plot_heatmap_gene_stage(dm_adata, gene_list: list[str]=None):
    """
    plot the expression of gene in gene_list over cell stages(E6.5 to E7.5).

    :param dm_adata: Anndata object
    :param gene_list: Genes to plot on the heatmap.
    :return: None(fig saved)
    """
    import numpy as np
    from anndata import AnnData
    from matplotlib.pyplot import show
    missing_genes = [gene for gene in gene_list if gene not in dm_adata.var_names]
    if missing_genes:
        raise ValueError(f"The following genes are not found in the dataset: {missing_genes}")

    adata_filtered = dm_adata[:, gene_list].copy()

    # 定义发育阶段到数值的映射
    stage_to_num = {"E65": 0, "E675": 1, "E70": 2, "E725": 3, "E75": 4}

    # 初始化存储平均表达量的矩阵
    n_genes = adata_filtered.shape[1]  # 基因数量
    n_stages = len(stage_to_num)  # 发育阶段数量
    mean_expression = np.zeros((n_stages, n_genes))  # 存储平均表达量

    # 计算每个发育阶段的平均表达量
    for stage, num in stage_to_num.items():
        # 提取当前发育阶段的细胞
        stage_cells = adata_filtered.obs["cell_time"] == stage
        # 计算当前发育阶段的基因平均表达量
        mean_expression[num, :] = np.mean(adata_filtered.X[stage_cells, :], axis=0)

    # 构建新的 AnnData 对象
    new_adata = AnnData(
        X=mean_expression,  # 平均表达量矩阵
        obs={"stage": list(stage_to_num.keys())},  # 发育阶段
        var=adata_filtered.var.copy(),  # 原有基因信息
    )

    # TODO：下面两个函数用着不方便，还是自己写一个heatmap函数吧
    sc.pl.timeseries_as_heatmap(new_adata.X, var_names=new_adata.var_names,
                                highlights_x=[], color_map="viridis")
    sc.pl._utils.savefig_or_show(writekey="heatmap_gene_stage", save=".png", show=False, dpi=300)


def setting():
    """
    do setting on if __name__ == "__main__": module in fa2_and_pseudotime.py

    :return: None
    """
    # 设置展示运行中会出现的信息
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    # 打印运行环境
    sc.logging.print_header()
    # 设置图像参数，保存信息
    sc.settings.set_figure_params(dpi=150, dpi_save=300, facecolor="white")
    sc.settings.figdir = SAVING_FIG_FOLDER
    sc.settings.autosave = IS_AUTOSAVE
    sc.settings.autoshow = IS_AUTOSHOW


if __name__ == "__main__":

    SAVING_FIG_FOLDER = './temp4'
    setting()
    adata = sc.read_h5ad("h5ads/final_627_combat_re.h5ad")
    cell_type_lst_a = ["Epiblast", "Primitive Streak", "Anterior Primitive Streak",
                       "Def. endoderm", "Gut", "Visceral endoderm", "ExE endoderm"]
    cell_type_lst_b = ["Epiblast", "Primitive Streak", "Nascent mesoderm", "Mixed mesoderm", "Mesenchyme",
                       "Haematoendothelial progenitors", "Blood progenitors 1", "Blood progenitors 2"]
    cell_type_lst_c = list(set(cell_type_lst_a + cell_type_lst_b))
    # NOTE: cell_50327 cell_79163 cell_86931 are outlier cells, filtered.
    # NOTE: cell_29635 cell_51875 cell_50514 cell_837, filtered.
    sub_adata = adata[adata.obs['cell_type'].isin(cell_type_lst_c), :]
    outlier_cells = ["cell_50327", "cell_79163", "cell_46746", "cell_29635", "cell_51875", "cell_50514", "cell_837",
                     "cell_86931"]
    sub_adata = sub_adata[~sub_adata.obs_names.isin(outlier_cells), :]
    sub_adata = paga(p_adata=sub_adata, color="cell_type", groups="cell_type")
    sub_adata = paga_scatter(ps_adata=sub_adata, color="cell_type")

    # sub_adata = pseudotime(sub_adata, n_branchings=1)
    with open("with_child_comp_genes.txt", "r") as file:
        gene_lst = [line.strip() for i, line in enumerate(file) if i < 400]
    # plot_pseudotime(sub_adata, gene_list=gene_lst, cell_sample_step=1)
    # up_gene_lst = ["Rbms1", "Car14", "Qprt", "Frzb", "Pmaip1", "Ccnd2"]
    # down_gene_lst = ["Cldn6", "Npm1", "Hnrnpa1", "Sgce", "Dap", "Stx7", "Ctsc", "Slc16a1", "Apoe", "Glrx", "Sms", "Nefl", "Slc25a4", "Ctsc"]
    # root_gene_lst = ["Dnmt3b", "Rab25", "Igfbp2", "Arl4c", "Smad1", "Bambi", "T", "Trh", "Lhx1", "Pdzd4"]
    prof_selected = ["Foxa2", "T", "Mesp1", "Lhx1", "Lhx6", "Gsc", "Hhex", "Cer1", "Phlda2", "Snai1", "Mixl1"]
    # plot_fa_color_pseudotime(sub_adata, color=up_gene_lst)
    # plot_fa_color_pseudotime(sub_adata, color=down_gene_lst)
    plot_fa_color_pseudotime(sub_adata, color=prof_selected)



    # 是否一定要先运行paga才再能运行paga_scatter呢？sc.tl.paga, sc.pl.paga, sc.tl.draw_graph, sc.pl.draw_graph函数逻辑复杂。
    # 实验证明不一定需要！因为在之前的绘制UMAP图时就跑过PAGA了。因为UMAP图的initial position就是PAGA的节点。sc.pl.draw_graph会自动判断传入的color参数是否为var中的基因，或对obs的注释。
    # TODO：
    #   如果为obs的categories类注释，则会自动细胞群体着色。如果为var_name中的基因，那么自动按照各单细胞中该基因的表达量着色。那么着色的值是什么？是scale后的那个表达量吗？
    #   而且draw_graph的得到的散点图的横轴纵轴的含义是？散点的分布和UMAP图的分布明显不一样，那么这又是什么算法？


