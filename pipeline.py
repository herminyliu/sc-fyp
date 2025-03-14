import scanpy as sc
import numpy as np
import preprocess
import umap_clustering
from anndata import AnnData

# 超参数定义
FILE_PATH = 'Expr_Mat_PS+Mesoderm_Label.txt'
# 要求细胞必须至少包含FILTER_CELLS_MIN_GENES种基因
preprocess.FILTER_CELLS_MIN_GENES = 200
# 要求基因必须至少出现在FILTER_GENES_MIN_CELLS个细胞中
preprocess.FILTER_GENES_MIN_CELLS = 5
# 只留下最多有MAX_GENES_BY_COUNTS种基因的单细胞，123代表批次
preprocess.MAX_GENES_BY_COUNTS_1 = 9700
preprocess.MAX_GENES_BY_COUNTS_2 = 11500
preprocess.MAX_GENES_BY_COUNTS_3 = 5700
PREPROCESS_SAVING_FILE_PATH = './preprocessed_data_123_noBatchEffect.h5ad'  # 应该为h5ad文件
CLUSTERING_SAVING_FILE_PATH = "./clustered_data_123_noBatchEffect.h5ad"
PREPROCESS_SAVING_FIG_FOLDER = './preprocessed_figures_123_noBatchEffect'
CLUSTERING_SAVING_FIG_FOLDER ='./clustering_figures_123_noBatchEffect'
# 建议IS_AUTOSAVE设置为False以避免图名相同覆盖问题
IS_AUTOSAVE = True
IS_AUTOSHOW = False
VERBOSITY = 2
N_MARKER_GENES = 20


def pipeline_settings(saving_fig_folder):
    # 设置展示运行中会出现的信息
    sc.settings.verbosity = VERBOSITY  # verbosity: errors (0), warnings (1), info (2), hints (3)
    # 打印运行环境
    sc.logging.print_header()
    # 设置图像参数，保存信息
    sc.settings.set_figure_params(dpi=150, dpi_save=300, facecolor="white")
    sc.settings.figdir = saving_fig_folder
    sc.settings.autosave = IS_AUTOSAVE
    sc.settings.autoshow = IS_AUTOSHOW


if __name__ == "__main__":
    # pipeline_settings(PREPROCESS_SAVING_FIG_FOLDER)
    # # 读取文本文件，并转置
    # adata = preprocess.read_origin_text(FILE_PATH)
    #
    # # 拆分行名，注释行
    # # adata = preprocess.annotate_cells(adata)
    #
    # adata.write_h5ad("./data_annotated_withBatch.h5ad")
    # # adata = sc.read_h5ad("./data_annotated_123.h5ad")
    #
    # # 检查数据中是否存在 NaN 值
    # if np.isnan(adata.X).any():
    #     print("Warning: Data contains NaN values. Handling them before proceeding.")
    # else:
    #     print("Data DO NOT contains NaN values.")
    # print("读取完成" + "=" * 50)
    # print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    #
    # # 检查是否有基因的同位体
    # preprocess.check_identical_cells(adata)
    #
    # # 数据预处理
    # adata = preprocess.preprocess_data(adata, with_batches=False)
    #
    # adata.write_h5ad(PREPROCESS_SAVING_FILE_PATH)
    # # adata = sc.read_h5ad(PREPROCESS_SAVING_FILE_PATH)
    #
    # print("预处理完成" + "=" * 50)
    # print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    # print(f"包含{adata.obs.cell_time.unique()} 时期的细胞")
    # print(f"包含{adata.obs.batch.unique()} batch的细胞")
    #
    pipeline_settings(CLUSTERING_SAVING_FIG_FOLDER)
    #
    # # 聚类分析
    # adata = umap_clustering.cluster_umap(adata, leiden_resolution=0.9)
    #
    # adata.write_h5ad(CLUSTERING_SAVING_FILE_PATH)

    adata = sc.read_h5ad(CLUSTERING_SAVING_FILE_PATH)

    # 绘制UMAP图像
    umap_clustering.plot_umap(adata)

    print("Clustering完成" + "=" * 50)
    print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    print(f"细胞被分成了{adata.obs.leiden.unique()}类")

    # 找Marker基因
    adata = umap_clustering.finding_marker_gene(adata, values_to_plot="logfoldchanges", groupby="cell_type", saving_fig_folder=CLUSTERING_SAVING_FIG_FOLDER, n_marker_gene=N_MARKER_GENES, method="t-test")
    adata = umap_clustering.finding_marker_gene(adata, values_to_plot="logfoldchanges", groupby="cell_time", saving_fig_folder=CLUSTERING_SAVING_FIG_FOLDER, n_marker_gene=N_MARKER_GENES, method="t-test")
    # adata = umap_clustering.finding_marker_gene(adata, groupby="leiden", saving_fig_folder=CLUSTERING_SAVING_FIG_FOLDER, n_marker_gene=N_MARKER_GENES, method="t-test")

    print("找marker基因完成" + "=" * 50)

    print("潜在marker基因dotplot绘图完成" + "=" * 50)

    # 依据Marker基因绘制针对基因的dotplot
    # adata = umap_clustering.plot_marker_gene(adata, marker_genes_lst=[])

    # 使用PAGA进行拟时序分析，PAGA
