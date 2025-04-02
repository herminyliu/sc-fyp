import os

import scanpy as sc
import numpy as np
# import mnnpy as mnn
from anndata import AnnData
from typing import Literal


# 超参数定义
FILE_PATH = './data_annotated_123.h5ad'
# 要求细胞必须至少包含FILTER_CELLS_MIN_GENES种基因
FILTER_CELLS_MIN_GENES = 200
# 要求基因必须至少出现在FILTER_GENES_MIN_CELLS个细胞中
FILTER_GENES_MIN_CELLS = 10
# 只留下最多有MAX_GENES_BY_COUNTS种基因的单细胞，123代表批次
MAX_GENES_BY_COUNTS_1 = 5000
MAX_GENES_BY_COUNTS_2 = 11020
MAX_GENES_BY_COUNTS_3 = 5337
PREPROCESS_SAVING_FILE_PATH = './preprocessed_data_123_noBatchEffect.h5ad'  # 应该为h5ad文件
PREPROCESS_SAVING_FIG_FOLDER = './preprocessed_figures_123_noBatchEffect'
IS_AUTOSAVE = False
IS_AUTOSHOW = True


def preprocess_data(pp_adata: AnnData, with_batches: bool=False, with_multistrip: bool=False,
                    method: Literal["combat", "mnn"]=None, is_regressout=True):
    """
    Preprocess single cell sequencing data in pp_adata.

    :param pp_adata: Single cell sequencing data, Anndata object.
    :param with_batches: Do the pp_adata contains obvious batch effect.
    :param with_multistrip: Do the scatter plot(x-axis is total_counts, all reads detected in one cell,
    yaxis is n_genes_by_counts, all types of genes detected in one cell) show more than one strips.
    :param method: The batch effect correction method to choose.
    :param is_regressout: Do scanpy.pp.regressout on the sequence depth or not.
    :return: Anndata object.
    """

    if method in ["combat", "mnn"] and with_batches is False:
        print("Warning: Claim no batch effect but still do batch correction. May lead to over smooth of the data!")

    # 绘制前n_top个最高表达基因的各单细胞内读数值分布箱线图
    sc.pl.highest_expr_genes(pp_adata, n_top=20)

    # 对细胞进行过滤，要求细胞必须至少包含min_genes种基因的reads
    sc.pp.filter_cells(pp_adata, min_genes=FILTER_CELLS_MIN_GENES, inplace=True)
    # 对细胞基因进行过滤，要求基因必须至少出现在min_cells个细胞中
    sc.pp.filter_genes(pp_adata, min_cells=FILTER_GENES_MIN_CELLS, inplace=True)

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
    if with_multistrip:
        sc.pl.scatter(pp_adata, x="total_counts", y="n_genes_by_counts", color="batch")
        # 只留下最多有MAX_GENES_BY_COUNTS种基因的单细胞
        bool_lst_1 = ((pp_adata.obs.batch == "1") & (pp_adata.obs.n_genes_by_counts < MAX_GENES_BY_COUNTS_1)) | (
                    (pp_adata.obs.batch == "2") & (pp_adata.obs.n_genes_by_counts < MAX_GENES_BY_COUNTS_2))
        bool_lst_2 = bool_lst_1 | (
                    (pp_adata.obs.batch == "3") & (pp_adata.obs.n_genes_by_counts < MAX_GENES_BY_COUNTS_3))
        pp_adata = pp_adata[bool_lst_2, :]
    else:
        sc.pl.scatter(pp_adata, x="total_counts", y="n_genes_by_counts")
        pp_adata = pp_adata[pp_adata.obs.n_genes_by_counts < MAX_GENES_BY_COUNTS_1, :]

    if pp_adata.n_obs == 0:
        raise ValueError("过滤后数据为空，请检查过滤条件！")

    # 将每个细胞的总读数归一化为target_sum指示的值
    sc.pp.normalize_total(pp_adata, target_sum=1e4, inplace=True, copy=False)

    # 将数据对数化
    sc.pp.log1p(pp_adata, copy=False)

    # 找到高变基因
    sc.pp.highly_variable_genes(pp_adata, n_top_genes=3000, flavor="seurat")

    # 对高变基因作图
    sc.pl.highly_variable_genes(pp_adata)

    if method == "combat":
        # DeepSeek建议“如果批次效应与某些技术性偏差（如测序深度）强相关，建议先去除这些技术性偏差，再处理批次效应”
        # Regress out (mostly) unwanted sources of variation.只对adata中的obs中的total_counts操作，方差为1
        # sc.pp.regress_out(pp_adata, keys=["batch"])
        if is_regressout:
            sc.pp.regress_out(pp_adata, keys=["total_counts"], copy=False)
        # 使用 ComBat 进行批次校正
        sc.pp.combat(pp_adata, key='batch', inplace=True)
    elif method == "mnn":
        # Do remember to receive the return value of sc.external.pp.mnn_correct!
        # It is recommended to pass log-transformed matrices/AnnData objects to mnn_correct, and use HVGs instead of all the genes.
        # If you use mnn, then no regress out will be performed.
        # The speed of mnn: Finishes correcting ~50000 cells/19 batches * ~30000 genes in ~12h on a 16 core 32GB mem server.
        if is_regressout:
            sc.pp.regress_out(pp_adata, keys=["total_counts"])
        pp_adata = pp_adata[:, pp_adata.var['highly_variable'] == True]
        adata_list = [pp_adata[pp_adata.obs['batch'] == batch].copy() for batch in pp_adata.obs['batch'].values.unique()]
        pp_adata, mnn_list, angle_list = mnn.mnn_correct(adata_list[0], adata_list[1], adata_list[2], batch_key='batch', var_subset=None,
                                   k=20, sigma=1.0, do_concatenate=True, n_jobs=16)
    elif not with_batches:
        if is_regressout:
            sc.pp.regress_out(pp_adata, keys=["total_counts"])
    else:
        print("Warning: Have not specified whether contain batch effect nor to use any method to correct!")
        print("Warning: Do as no batch effect and do no correction!")
        if is_regressout:
            sc.pp.regress_out(pp_adata, keys=["total_counts"])

    if not isinstance(pp_adata, AnnData):
        raise ValueError(f"In preprocess.py, after batch correction, pp_adata no longer is an AnnData object. The type is{type(pp_adata)}. The value is{pp_adata}")

    # 将数据归一化，均值为0，方差为1，归一化后大于max_value的值截断到max_value, zero_center这个参数代表是否要将数据的均值归到0
    # 如果zero_center为True，那么数据会只保留[-max_value, max_value]
    # 如果zero_center为False，那么数据会只保留[-inf, max_value]
    sc.pp.scale(pp_adata, max_value=10, zero_center=True)

    return pp_adata


def annotate_cells(a_adata: AnnData):
    """
    Do cell annotation based on the row(observation) index.

    :param a_adata: Single cell sequencing data, Anndata object.
    :return: Anndata object.
    """
    # 这里必须要注意，不能写成adata_1 = adata[adata.obs.index.str.split('.')[-2] == 1, :]
    # adata.obs.index.str.split('.')返回的是一个“列表的列表”
    # .str[-2]代表的操作是“对列表中的每一个元素（这是列表），取每个元素（这是列表）中的倒数第二个元素，并且替代（注意替代）原有元素”
    # 首先，需要和字符串1作比较；其次，最后一个.str是pandas提供的针对series和dataframe的字符串方法，可以逐元素操作。

    # Retrieve batch info.
    batch_info = a_adata.obs.index.str.split('.').str[-2]

    # 提取老师预先注释好的细胞类别信息
    celltype_info = a_adata.obs.index.str.split('.').str[-1]

    # 提取细胞发育时段信息
    celltime_info = a_adata.obs.index.str.split('.').str[-3]

    # 提取细胞的ID
    cellID_info = a_adata.obs.index.str.split('.').str[0].str.split('_').str[1]

    # 提取细胞的胚胎ID
    embryo_info = a_adata.obs.index.str.split('.').str[1]

    # 将信息添加到 obs 中
    a_adata.obs['batch'] = batch_info.astype('category')
    a_adata.obs['cell_type'] = celltype_info.astype('category')
    a_adata.obs['cell_time'] = celltime_info.astype('category')
    a_adata.obs['cell_ID'] = cellID_info
    a_adata.obs['embryo_ID'] = embryo_info.astype('category')
    a_adata.obs["cell_type_abbr"] = a_adata.obs["cell_type"].cat.rename_categories(
        {
            "Anterior_Primitive_Streak": "APS",
            "Intermediate_mesoderm": "IM",
            "Mixed_mesoderm": "MM",
            "Nascent_mesoderm": "NM",
            "Paraxial_mesoderm": "PAM",
            "Pharyngeal_mesoderm": "PHM",
            "Primitive_Streak": "PS",
            "Somitic_mesoderm": "SM",
        }
    )
    return a_adata


def annotate_genes(a_adata: AnnData, gene_id_mapping_file_path: str):
    """
    Replace gene IDs in the AnnData object's var (columns) with gene names based on a TSV file.

    Parameters:
    - a_adata: AnnData object where the var (columns) contain gene IDs.
    - gene_id_file_path: str, path to the TSV file containing gene IDs and their corresponding gene names.

    The TSV file should have the following format:
    - Two columns: the first column is the gene ID, and the second column is the gene name.
    - No header row.
    - Example:
        ENSMUSG00000051951    Xkr4
        ENSMUSG00000051952    Rp1
        ENSMUSG00000051953    Sox17

    Returns:
    - Anndata object with annotated and deduplicated genes.
    """
    from pandas import read_csv
    import anndata as ad
    # Read the TSV file into a DataFrame
    gene_id_df = read_csv(gene_id_mapping_file_path, sep='\t', header=None, names=['gene_id', 'gene_name'])

    # Create a dictionary mapping gene IDs to gene names
    gene_id_to_name = dict(zip(gene_id_df['gene_id'], gene_id_df['gene_name']))

    # Replace gene IDs in the AnnData object's var with gene names
    # Only replace gene IDs that exist in the mapping dictionary
    a_adata.var.index = a_adata.var.index.map(lambda x: gene_id_to_name.get(x, x))

    # Print a summary of the replacement
    replaced_count = sum(a_adata.var.index.isin(gene_id_to_name.values()))
    print(f"Replaced {replaced_count} gene IDs with gene names. The total number of genes is {a_adata.n_vars}")
    # NOTE: Very important! Do remember to check duplicates after replace the IDs into the actual gene names.
    if a_adata.var.index.is_unique:
        print("adata.var.index do not contain duplicated genes.")
        return a_adata
    else:
        # Identify duplicate indices in var.index
        duplicate_gene_names = a_adata.var.index[a_adata.var.index.duplicated(keep='first')]

        # Create a dictionary to store the sum of duplicate columns
        sum_dict = {}

        # Iterate over duplicate indices and sum the corresponding columns
        for gene_name in duplicate_gene_names:
            # Select columns with the same index
            duplicate_columns = a_adata.X[:, a_adata.var.index == gene_name]
            # Sum the columns along the axis (axis=1 for rows, axis=0 for columns)
            sum_dict[gene_name] = np.sum(duplicate_columns, axis=1)

        # Keep a copy of the duplicated
        # To keep all the complicated obs and var annotations
        a_adata_duplicated = a_adata[:, a_adata.var.index.duplicated(keep="first")]

        # NOTE: There seems a bug in pandas method pandas.Index.duplicated.
        # NOTE: a_adata.var.index.duplicated(keep="first") should only keep the index of duplicated genes first shown.
        # NOTE: But if a gene has three duplication, this method will keep two copies instead of one!

        # if not a_adata_duplicated.var.index.is_unique:
        #     print(f"Warning: a_adata_duplicated.var.index contain duplicated gene names."
        #                      f"a_adata_duplicated.var.index:{a_adata_duplicated.var.index}")
        #     a_adata_duplicated = a_adata_duplicated[:, ~a_adata_duplicated.var.index.duplicated(keep="first")]

        def remove_duplicated_genes(duplicated: AnnData) -> AnnData:
            """
            Recursively remove duplicate gene names in a_adata_duplicated.var.index, until no duplicates

            :param duplicated: anndata.AnnData object, containing gene expression data.

            :return: anndata.AnnData: deduplicated AnnData object，making sure no duplicates in var.index.
            """
            if not duplicated.var.index.is_unique:
                print(f"Warning: a_adata_duplicated.var.index contains duplicated gene names. "
                      f"a_adata_duplicated.var.index: {duplicated.var.index}")
                # 移除重复值，保留第一个出现的基因名
                duplicated = duplicated[:, ~duplicated.var.index.duplicated(keep="first")]
                # 递归调用，直到没有重复值
                return remove_duplicated_genes(duplicated)
            else:
                # 如果没有重复值，返回处理后的 AnnData 对象
                return duplicated

        a_adata_duplicated = remove_duplicated_genes(a_adata_duplicated)

        # Remove ALL duplicate columns from the AnnData object
        a_adata = a_adata[:, ~a_adata.var.index.duplicated(keep=False)]

        for gene_name, summed_vector in sum_dict.items():
            a_adata_duplicated.X[:, a_adata_duplicated.var.index.get_loc(gene_name)] = summed_vector

        b_adata = ad.concat(adatas=[a_adata, a_adata_duplicated], axis='var', join='outer', merge="same", uns_merge=None)

        if b_adata.var_names is None or b_adata.obs_names is None:
            raise ValueError("b_adata.var_names is None or b_adata.obs_names is None")
        return b_adata


def annotate_cell_id(a_adata: AnnData):
    """
    annotate ID of cells

    :param a_adata: anndata.AnnData
    :return: anndata.AnnData
    """
    cellID_info = a_adata.obs.index.str.split('_').str[1]
    a_adata.obs['cell_ID'] = cellID_info
    return a_adata


def transfer_annotations(src_adata: AnnData, dest_adata: AnnData):
    """
    Transfer the annotations from src_adata to dest_adata one-to-one based on cell_ID.

    Params:
        src_adata: reference AnnData object.
        dest_adata: target AnnData object, accepting annotation from reference.
    """
    # 确保 cell_ID 存在于两个 AnnData 对象中
    if "cell_ID" not in src_adata.obs.columns or "cell_ID" not in dest_adata.obs.columns:
        raise ValueError("Both src_adata and dest_adata must have 'cell_ID' in obs.")

    # 获取 src_adata 中的注释列
    annotations = ["batch", "cell_type", "cell_time", "embryo_ID", "cell_type_abbr"]
    for col in annotations:
        if col not in src_adata.obs.columns:
            raise ValueError(f"Column '{col}' not found in src_adata.obs.")

    # 将 src_adata.obs 转换为字典，以 cell_ID 为键
    src_obs_dict = src_adata.obs.set_index("cell_ID").to_dict(orient="index")

    # 遍历 dest_adata 的 cell_ID，将注释从 src_adata 转移到 dest_adata
    for cell_id in dest_adata.obs["cell_ID"]:
        if cell_id in src_obs_dict:
            for col in annotations:
                dest_adata.obs.loc[dest_adata.obs["cell_ID"] == cell_id, col] = src_obs_dict[cell_id][col]
        else:
            # 如果 cell_ID 在 src_adata 中不存在，填充为 NaN
            for col in annotations:
                dest_adata.obs.loc[dest_adata.obs["cell_ID"] == cell_id, col] = None

    print("Annotations transferred successfully!")

    return dest_adata


def check_identical_cells(u_adata: AnnData):
    """
    check whether u_adata contains duplicated rows(cells).

    :param u_adata: anndata.AnnData
    :return: anndata.AnnData
    """
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
    """
    setting hyper params, only effective when directly run preprocess.py. This function only used in
    if __name__ == "__main__" module of preprocess.py

    :return: None
    """
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
    """
    Read original data in txt format.
    Note: since the txt file takes cells as columns while genes as rows, a transpose action is applied. If your dataset
    Note: do not follow this rule, do remove transposed_X = origin_adata.X.T

    :param file_path: txt file path
    :return: anndata.AnnData
    """
    origin_adata = sc.read_text(file_path, delimiter='\t')
    transposed_X = origin_adata.X.T
    adata = sc.AnnData(X=transposed_X, var=origin_adata.obs, obs=origin_adata.var)
    return adata


if __name__ == '__main__':
    # Test whether preprocess.py run normally.
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
