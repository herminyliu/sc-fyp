import scanpy as sc
import numpy as np
from anndata import AnnData

# 超参数定义
FILE_PATH = './data_annotated_123.h5ad'
# 要求细胞必须至少包含FILTER_CELLS_MIN_GENES种基因
FILTER_CELLS_MIN_GENES = 200
# 要求基因必须至少出现在FILTER_GENES_MIN_CELLS个细胞中
FILTER_GENES_MIN_CELLS = 10
# 只留下最多有MAX_GENES_BY_COUNTS种基因的单细胞，123代表批次
MAX_GENES_BY_COUNTS_1 = 9700
MAX_GENES_BY_COUNTS_2 = 11500
MAX_GENES_BY_COUNTS_3 = 5700
PREPROCESS_SAVING_FILE_PATH = './preprocessed_data_123_noBatchEffect.h5ad'  # 应该为h5ad文件
PREPROCESS_SAVING_FIG_FOLDER = './preprocessed_figures_123_noBatchEffect'
IS_AUTOSAVE = False
IS_AUTOSHOW = True


def preprocess_data(pp_adata, with_batches=False):
    # pp代表preprocessing
    # 对细胞进行细胞发育时期等注释
    pp_adata = annotate_cells(pp_adata)

    # 绘制前n_top个最高表达基因的各单细胞内读数值分布箱线图
    sc.pl.highest_expr_genes(pp_adata, n_top=20)

    # 对细胞进行过滤，要求细胞必须至少包含min_genes种基因的reads
    sc.pp.filter_cells(pp_adata, min_genes=FILTER_CELLS_MIN_GENES)
    # 对细胞基因进行过滤，要求基因必须至少出现在min_cells个细胞中
    sc.pp.filter_genes(pp_adata, min_cells=FILTER_GENES_MIN_CELLS)

    # 计算质量控制指标，inplace=True代表函数会添加varm/obsm/varp/obsp（细胞注释或基因注释）到adata中
    sc.pp.calculate_qc_metrics(
        pp_adata, percent_top=None, log1p=False, inplace=True
    )

    # 绘制小提琴图，一共两幅
    # 第一幅n_genes_by_counts，绘制单细胞内基因种类数目分布情况
    # 第二幅total_counts，绘制所有单细胞的读数的值分布情况
    sc.pl.violin(
        pp_adata,
        keys= ["n_genes_by_counts", "total_counts"],
        jitter=0.4,
        multi_panel=True,
    )

    # 绘制散点图，横纵轴含义如x，y所示，均为adata.obs
    if with_batches:
        sc.pl.scatter(pp_adata, x="total_counts", y="n_genes_by_counts", color="batch")
    else:
        sc.pl.scatter(pp_adata, x="total_counts", y="n_genes_by_counts")

    # 只留下最多有MAX_GENES_BY_COUNTS种基因的单细胞
    bool_lst_1 = ((pp_adata.obs.batch == "1") & (pp_adata.obs.n_genes_by_counts < MAX_GENES_BY_COUNTS_1)) | ((pp_adata.obs.batch == "2") & (pp_adata.obs.n_genes_by_counts < MAX_GENES_BY_COUNTS_2))
    bool_lst_2 = bool_lst_1 | ((pp_adata.obs.batch == "3") & (pp_adata.obs.n_genes_by_counts < MAX_GENES_BY_COUNTS_3))

    pp_adata = pp_adata[bool_lst_2, :]
    if pp_adata.n_obs == 0:
        raise ValueError("过滤后数据为空，请检查过滤条件！")

    # 将每个细胞的总读数归一化为target_sum指示的值
    sc.pp.normalize_total(pp_adata, target_sum=1e4)

    # 将数据对数化
    sc.pp.log1p(pp_adata)

    if with_batches:
        # DeepSeek建议“如果批次效应与某些技术性偏差（如测序深度）强相关，建议先去除这些技术性偏差，再处理批次效应”
        # Regress out (mostly) unwanted sources of variation.只对adata中的obs中的total_counts操作，方差为1
        sc.pp.regress_out(pp_adata, keys=["batch"])
        sc.pp.regress_out(pp_adata, keys=["total_counts"])
        # 使用 ComBat 进行批次校正
        sc.pp.combat(pp_adata, key='batch')
    else:
        sc.pp.regress_out(pp_adata, keys=["total_counts"])

    # 找到高变基因
    sc.pp.highly_variable_genes(pp_adata, n_top_genes=3000, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor="seurat")

    # 对高变基因作图
    sc.pl.highly_variable_genes(pp_adata)

    # 将数据归一化，均值为0，方差为1，归一化后大于max_value的值截断到max_value, zero_center这个参数代表是否要将数据的均值归到0
    # 如果zero_center为True，那么数据会只保留[-max_value, max_value]
    # 如果zero_center为False，那么数据会只保留[-inf, max_value]
    sc.pp.scale(pp_adata, max_value=10, zero_center=True)

    return pp_adata


def annotate_cells(a_adata: AnnData):
    # 这里必须要注意，不能写成adata_1 = adata[adata.obs.index.str.split('.')[-2] == 1, :]
    # adata.obs.index.str.split('.')返回的是一个“列表的列表”
    # .str[-2]代表的操作是“对列表中的每一个元素（这是列表），取每个元素（这是列表）中的倒数第二个元素，并且替代（注意替代）原有元素”
    # 首先，需要和字符串1作比较；其次，最后一个.str是pandas提供的针对series和dataframe的字符串方法，可以逐元素操作。

    # 提取批次信息
    batch_info = a_adata.obs.index.str.split('.').str[-2]

    # 提取老师预先注释好的细胞类别信息
    celltype_info = a_adata.obs.index.str.split('.').str[-1]

    # 提取细胞发育时段信息
    celltime_info = a_adata.obs.index.str.split('.').str[-3]

    # 提取细胞的ID
    cellID_info = a_adata.obs.index.str.split('.').str[0]

    # 将信息添加到 obs 中
    a_adata.obs['batch'] = batch_info.astype('category')
    a_adata.obs['cell_type'] = celltype_info.astype('category')
    a_adata.obs['cell_time'] = celltime_info.astype('category')
    a_adata.obs['cell_ID'] = cellID_info
    return a_adata


def check_identical_cells(u_adata: AnnData):
    from scipy.sparse import csr_matrix
    import numpy as np

    # 检查数据是否是稀疏矩阵
    if isinstance(u_adata.X, csr_matrix):
        data = u_adata.X.toarray()
        print("储存格式为稀疏矩阵")
    else:
        data = u_adata.X
        print("储存格式不为稀疏矩阵")

    # 检查是否存在完全相同的细胞
    unique_cells, indices = np.unique(data, axis=0, return_inverse=True)
    if len(unique_cells) < len(data):
        print(f'Found {len(data) - len(unique_cells)} duplicate cells.')
    else:
        print('No duplicate cells found.')

    # u_adata = sc.pp.scrublet(u_adata, copy=True, batch_key='batch')
    # print("The number of predicted scrublet cells is:", np.sum(u_adata.obs["predicted_doublet"] == True))


def setting():
    # 设置展示运行中会出现的信息
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    # 打印运行环境
    sc.logging.print_header()
    # 设置图像参数，保存信息
    sc.settings.set_figure_params(dpi=150, dpi_save=300, facecolor="white")
    sc.settings.figdir = PREPROCESS_SAVING_FIG_FOLDER
    sc.settings.autosave = IS_AUTOSAVE
    sc.settings.autoshow = IS_AUTOSHOW


def read_origin_text(file_path):
    origin_adata = sc.read_text(file_path, delimiter='\t')
    transposed_X = origin_adata.X.T
    adata = sc.AnnData(X=transposed_X, var=origin_adata.obs, obs=origin_adata.var)
    return adata


if __name__ == '__main__':
    setting()
    # 读取文本文件，并转置
    adata = read_origin_text(PREPROCESS_SAVING_FILE_PATH)

    adata = annotate_cells(adata)

    adata.write_h5ad("./data_annotated_123.h5ad")

    # 检查数据中是否存在 NaN 值
    if np.isnan(adata.X).any():
        print("Warning: Data contains NaN values. Handling them before proceeding.")
    else:
        print("Data DO NOT contains NaN values.")

    check_identical_cells(adata)

    adata = preprocess_data(adata, with_batches=True)

    print("预处理完成" + "="*50)
    print(f"细胞数: {adata.n_obs}, 基因数: {adata.n_vars}")
    print(f"包含{adata.obs.cell_time.unique()} 时期的细胞")
    print(f"包含{adata.obs.batch.unique()} batch的细胞")
