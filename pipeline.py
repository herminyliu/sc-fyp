import scanpy as sc
import numpy as np
import preprocess
import select_gene
import fa2_and_pseudotime
import umap_clustering


TXT_PATH_1 = './raw_txt/Expr_Mat_PS+Mesoderm_Filter_ComBat_No.Rep_Cell.Name.Merge_2.txt'
TXT_PATH_2 = './raw_txt/Expr_Mat_PS+Mesoderm_Label.txt'
# 要求细胞必须至少包含FILTER_CELLS_MIN_GENES种基因
preprocess.FILTER_CELLS_MIN_GENES = 500
# 要求基因必须至少出现在FILTER_GENES_MIN_CELLS个细胞中
preprocess.FILTER_GENES_MIN_CELLS = 5
# 只留下最多有MAX_GENES_BY_COUNTS种基因的单细胞，123代表批次
# Batch 1: 前 5% 最大值 = 9534.0
# Batch 1: 前 10% 最大值 = 9477.0
# Batch 2: 前 5% 最大值 = 11020.0
# Batch 2: 前 10% 最大值 = 10988.0
# Batch 3: 前 5% 最大值 = 5337.0
# Batch 3: 前 10% 最大值 = 5157.0
preprocess.MAX_GENES_BY_COUNTS_1 = 5000
preprocess.MAX_GENES_BY_COUNTS_2 = 11020
preprocess.MAX_GENES_BY_COUNTS_3 = 5337

SAVE_FIG_FOLDER = './to_blood_627_combat_figures_re'
FINAL_SAVE_H5AD_FILE = "./h5ads/final_to_blood_627_combat_re.h5ad"
PREPROCESSED_SAVE_H5AD_FILE = "./h5ads/raw_627_combat.h5ad"
# 建议IS_AUTOSAVE设置为False以避免图名相同覆盖问题
IS_AUTOSAVE = True
IS_AUTOSHOW = False
VERBOSITY = 2
N_MARKER_GENES = 20


def pipeline_settings(saving_fig_folder):
    """
    Set the hyper params in the python script pipeline.py

    :param saving_fig_folder: The folder where the store the figures drawn in this script.
    NOTE: this param can only control the figures drawn by scanpy.pl
    :return: No return
    """
    # 设置展示运行中会出现的信息
    sc.settings.verbosity = VERBOSITY  # verbosity: errors (0), warnings (1), info (2), hints (3)
    # 打印运行环境
    sc.logging.print_header()
    # 设置图像参数，保存信息
    sc.settings.set_figure_params(dpi=150, dpi_save=300, facecolor="white")
    sc.settings.figdir = saving_fig_folder
    sc.settings.autosave = IS_AUTOSAVE
    sc.settings.autoshow = IS_AUTOSHOW


def check_data(c_adata):
    """
    Check the dimension and annotation of c_adata.

    :param c_adata: Anndata object
    :return: No return
    """
    print(f"细胞数: {c_adata.n_obs}, 基因数: {c_adata.n_vars}")
    print(f"细胞类别有：{c_adata.obs.cell_type.unique()}")
    print(f"对细胞的注释有:{c_adata.obs.columns.unique()}")
    print(f"对基因的注释有:{c_adata.var.columns.unique()}")
    print(f"对细胞的多维注释有:{c_adata.obsm_keys()}")
    print(f"对基因的多维注释有:{c_adata.varm_keys()}")
    if len(c_adata.obs.cell_time.unique()) < 5:
        raise ValueError(f"只包含了{c_adata.obs.cell_time.unique()} 时期的细胞")
    if len(c_adata.obs.batch.unique()) < 3:
        raise ValueError(f"只包含了{c_adata.obs.batch.unique()} batch的细胞")


if __name__ == "__main__":
    pipeline_settings(SAVE_FIG_FOLDER)
    # # 读取文本文件，并转置
    # src_adata = preprocess.read_origin_text(TXT_PATH_1)
    # dest_adata = preprocess.read_origin_text(TXT_PATH_2)
    #
    # # 拆分行名，注释行
    # src_adata = preprocess.annotate_cells(src_adata)
    # dest_adata = preprocess.annotate_cell_id(dest_adata)
    #
    # adata = preprocess.transfer_annotations(src_adata=src_adata, dest_adata=dest_adata)
    # check_data(adata)
    #
    # # --------------------------
    adata = sc.read_h5ad(PREPROCESSED_SAVE_H5AD_FILE)
    adata.obs.rename(columns={"stage": "cell_time", "cell": "cell_ID", "sequencing.batch": "batch", "celltype": "cell_type", "sample": "embryo_ID"}, inplace=True)

    adata = preprocess.annotate_genes(adata, gene_id_mapping_file_path=r"E:\Studying\atlas\genes.tsv")
    #
    print("=" * 50+"对基因的注释完成"+"=" * 50)
    # check_data(adata)
    #
    # # 18号胚胎可能存在明显批次效应，删除
    # adata = adata[adata.obs["embryo_ID"] != str(18)]
    #
    # # 选取部分细胞类型
    # cell_type_lst_a = ["Epiblast", "Primitive Streak", "Anterior Primitive Streak",
    #                    "Def. endoderm", "Gut", "Visceral endoderm", "ExE endoderm"]
    # outlier_cells = ["cell_50327", "cell_79163", "cell_46746", "cell_29635", "cell_51875", "cell_50514", "cell_837",
    #                  "cell_86931"]
    # adata = adata[adata.obs['cell_type'].isin(cell_type_lst_a), :]
    #
    # check_data(adata)
    # adata = adata[~adata.obs_names.isin(outlier_cells), :]
    #
    #
    # # 检查数据中是否存在 NaN 值
    # # if np.isnan(adata.X).any():
    # #     print("Warning: Data contains NaN values. Handling them before proceeding.")
    # # else:
    # #     print("Data DO NOT contains NaN values.")
    # # print("=" * 50+"读取完成" + "=" * 50)
    # check_data(adata)
    #
    # # 检查是否有基因的同位体
    # preprocess.check_identical_cells(u_adata=adata)
    #
    # # 数据预处理
    # adata = preprocess.preprocess_data(adata, with_batches=False, method="combat", with_multistrip=False, is_regressout=True)
    #
    # adata.write_h5ad(PREPROCESSED_SAVE_H5AD_FILE)

    # 选取部分细胞类型
    cell_type_lst_a = ["Epiblast", "Primitive Streak", "Anterior Primitive Streak",
                       "Def. endoderm", "Gut", "Visceral endoderm", "ExE endoderm"]
    cell_type_lst_b = ["Epiblast", "Primitive Streak", "Nascent mesoderm", "Mixed mesoderm", "Mesenchyme",
                       "Haematoendothelial progenitors", "Blood progenitors 1", "Blood progenitors 2"]
    cell_type_lst_c = list(set(cell_type_lst_a + cell_type_lst_b))

    outlier_cells = ["cell_50327", "cell_79163", "cell_46746", "cell_29635", "cell_51875", "cell_50514", "cell_837",
                     "cell_86931"]
    adata = adata[adata.obs['cell_type'].isin(cell_type_lst_c), :]
    adata = adata[~adata.obs_names.isin(outlier_cells), :]

    print("=" * 50+"预处理完成" + "=" * 50)
    check_data(adata)

    # 聚类分析
    adata = umap_clustering.cluster_umap(adata, leiden_resolution=0.9, n_neighbors=15, spread=1.0, min_dist=0.5)
    print("=" * 50+"UMAP聚类完成" + "=" * 50)
    # 绘制UMAP图像
    umap_clustering.plot_umap(adata)
    print("=" * 50+"UMAP图绘制完成"+"=" * 50)
    check_data(adata)
    print(f"细胞被分成了{adata.obs.leiden.unique()}类")

    # 找Marker基因
    umap_clustering.finding_marker_gene(adata, values_to_plot=None, groupby="cell_type", saving_fig_folder=SAVE_FIG_FOLDER, n_marker_gene=N_MARKER_GENES, method="t-test")
    umap_clustering.finding_marker_gene(adata, values_to_plot=None, groupby="cell_time", saving_fig_folder=SAVE_FIG_FOLDER, n_marker_gene=N_MARKER_GENES, method="t-test")
    umap_clustering.finding_marker_gene(adata, values_to_plot=None, groupby="leiden", saving_fig_folder=SAVE_FIG_FOLDER, n_marker_gene=N_MARKER_GENES, method="t-test")

    print("=" * 50+"marker基因dotplot图绘制完成"+"=" * 50)

    # New gene of interest list mnn.
    gene_lst = ["Lhx1", "Lefty2", "Cyp26a1", "Fgf3", "Phlda2", "Bex4", "Epcam", "Trh", "Dnmt3b", "Pim2"]
    # 依据Marker基因绘制针对基因的dotplot
    umap_clustering.plot_marker_gene(adata, marker_genes_lst=gene_lst)
    check_data(adata)

    # 使用PAGA进行拟时序分析
    adata = fa2_and_pseudotime.paga(p_adata=adata, color=['embryo_ID', "cell_type", "cell_time"], groups="cell_type")
    adata = fa2_and_pseudotime.paga(p_adata=adata, color=['embryo_ID', "cell_type", "cell_time"], groups="leiden")
    adata = fa2_and_pseudotime.paga(p_adata=adata, color=['embryo_ID', "cell_type", "cell_time"], groups="cell_time")
    print("=" * 50+"PAGA拟时序分析完成，PAGA抽象图保存" + "=" * 50)
    adata = fa2_and_pseudotime.paga_scatter(ps_adata=adata, color=['embryo_ID', "cell_type", "cell_time"])
    print("=" * 50+"PAGA拟时序分析完成，ForceAtlas2图保存" + "=" * 50)
    check_data(adata)

    # 计算扩散图
    adata = fa2_and_pseudotime.diff(adata, color=['embryo_ID', "cell_type", "cell_time"])
    print("=" * 50+ "扩散分析完成" + "=" * 50)

    # 计算拟时间
    adata = fa2_and_pseudotime.pseudotime(dm_adata=adata, n_branchings=1)
    # gene_lst = select_gene.select_genes_by_dpt_trend(adata=adata, n_genes=10)
    print("=" * 50 + "DPT拟时序分析完成"+"=" * 50)
    check_data(adata)
    adata = sc.read_h5ad("h5ads/final_627_combat_re.h5ad")

    with open("with_child_comp_genes.txt", "r") as file:
        gene_lst = [line.strip() for i, line in enumerate(file) if i < 20]

    # 绘图
    fa2_and_pseudotime.plot_pseudotime(dm_adata=adata, gene_list=gene_lst, cell_sample_step=5)
    print("=" * 50+ "DPT拟时序基因图绘制完成" + "=" * 50)
    fa2_and_pseudotime.plot_fa_color_pseudotime(adata)
    print("=" * 50+ "按DPT拟时序着色的ForceAtlas图绘制完成" + "=" * 50)

    fa2_and_pseudotime.plot_heatmap_gene_stage(adata, gene_list=gene_lst)
    print("=" * 50+ "热图完成" + "=" * 50)
    adata.write_h5ad(FINAL_SAVE_H5AD_FILE)
