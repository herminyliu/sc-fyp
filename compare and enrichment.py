import os

import scanpy as sc
import pandas as pd

import fa2_and_pseudotime


def find_transition_genes(
        adata: sc.AnnData,
        cell_group_a: str,
        cell_group_b: str,
        groupby: str = "cell_type",
        min_logfc: float = 0.25,
        max_padj: float = 0.05,
        key_added: str = "transition_genes",
        save_csv: bool = True
) -> pd.DataFrame:
    """
    识别驱动细胞从类型A向类型B转变的差异表达基因

    参数：
        adata: 经过预处理的AnnData对象
        cell_type_a: 起始细胞类型
        cell_type_b: 目标细胞类型
        groupby: obs列中存储细胞类型信息的列名
        min_logfc: 最小对数倍数变化阈值
        max_padj: 最大校正后p值
        key_added: 结果存储在uns中的键名
        save_csv: 是否保存csv文件，会保存两个文件，分别记录基因名和对应的p值，需要将两个文件合并起来看。

    返回：
        包含显著差异基因及统计量的DataFrame
    """
    # 校验输入
    if groupby not in adata.obs:
        raise ValueError(f"Column '{groupby}' not found in adata.obs_names:{adata.obs}")

    cell_types = adata.obs[groupby].unique()
    if cell_group_a not in cell_types:
        raise ValueError(f"Cell type '{cell_group_a}' not found")
    if cell_group_b not in cell_types:
        raise ValueError(f"Cell type '{cell_group_b}' not found")

    # 创建细胞子集
    mask = adata.obs[groupby].isin([cell_group_a, cell_group_b])
    adata_sub = adata[mask].copy()

    # 设置比较组
    adata_sub.obs["comparison_group"] = [
        "target" if ct == cell_group_b else "reference"
        for ct in adata_sub.obs[groupby]
    ]

    # 进行差异表达分析
    sc.tl.rank_genes_groups(
        adata_sub,
        groupby="comparison_group",
        groups=["target"],  # params group accept list object.
        reference="reference",
        method="wilcoxon",
        use_raw=False  # 使用已经归一化后的数据
    )

    from pandas import DataFrame
    if save_csv:
        DataFrame(adata_sub.uns['rank_genes_groups']['pvals']).to_csv(f"./csv/discrimination_{groupby}_{cell_group_a}_{cell_group_b}_pvals.csv")
        DataFrame(adata_sub.uns['rank_genes_groups']['names']).to_csv(f"./csv/discrimination_{groupby}_{cell_group_a}_{cell_group_b}_names.csv")
        print(f"./csv/discrimination_{groupby}_{cell_group_a}_{cell_group_b}_pvals.csv successful saved.")


    # 提取结果
    result = sc.get.rank_genes_groups_df(adata_sub, group="target")

    # 过滤显著基因
    significant_genes = result[
        ((result["logfoldchanges"] >= min_logfc) | (result["logfoldchanges"] <= -min_logfc)) &
        (result["pvals_adj"] <= max_padj)
        ].sort_values("logfoldchanges", ascending=False)
    result = result.sort_values(by="logfoldchanges", ascending=False)

    # 存储结果到原始adata对象
    adata.uns[key_added] = significant_genes

    # 保存显著基因df
    if save_csv:
        significant_genes.to_csv(f"./csv/discrimination_{groupby}_{cell_group_a}_{cell_group_b}_significant_genes.csv")
        print(f"./csv/discrimination_{groupby}_{cell_group_a}_{cell_group_b}_significant_genes.csv successful saved.")

    return result


def plot_volcano(result_df,
                 logfc_thresh=0.5,
                 padj_thresh=0.05,
                 cell_group_a="Primitive Streak",
                 cell_group_b="Nascent mesoderm",
                 save_dir="discrimination_plot"):
    """
    Draw Volcano Plot.Y-axis is -Log10(Adjusted p-value), X-axis is logfoldchange. Downregulated genes drawn on the left half in blue, while the upregulated genes drawn on
    the right half in red.

    :param result_df: pandas.DataFrame 其中含有cell_group_a, cell_group_b的基因数据.
    每一行为基因，有多列，在本代码中会使用'logfoldchanges', 'pvals_adj'等列
    :param logfc_thresh: The minium logfoldchange value to accept, bigger than this value will be filtered.
    :param padj_thresh: The maximum p-value to accept, bigger than this value will be filtered.
    :param cell_group_a: cell_group_a: cell group a grouped by the criteria groupby
    :param cell_group_b: cell_group_b: cell group a grouped by the criteria groupby
    :param save_dir: The dictionary where the figs to be stored, ALL csv files in this repo go directly to ./csv
    :return: No return. Just save the fig.
    """
    import matplotlib.pyplot as plt
    from numpy import log10
    plt.figure(figsize=(8, 6))

    # 创建颜色分类列significance，用于散点图着色
    result_df['significance'] = 'n.s.'
    result_df.loc[(result_df['logfoldchanges'] >= logfc_thresh) &
                  (result_df['pvals_adj'] <= padj_thresh), 'significance'] = 'Up'
    result_df.loc[(result_df['logfoldchanges'] <= -logfc_thresh) &
                  (result_df['pvals_adj'] <= padj_thresh), 'significance'] = 'Down'

    # result_df['pvals_adj']中有部分基因的值非常低，直接记为0，会导致对数转换错误
    result_df['pvals_adj'].apply(lambda x: 1 * 10 ** (-300) if x == 0 else x)

    # 绘制散点
    scatter = plt.scatter(
        x=result_df['logfoldchanges'],
        y=-log10(result_df['pvals_adj']),
        c=result_df['significance'].map({'n.s.': 'grey', 'Up': 'red', 'Down': 'blue'}),
        s=10,
        alpha=0.6
    )

    # 添加阈值线
    plt.axvline(logfc_thresh, color='black', linestyle='--', linewidth=0.8)
    plt.axvline(-logfc_thresh, color='black', linestyle='--', linewidth=0.8)
    plt.axhline(-log10(padj_thresh), color='black', linestyle='--', linewidth=0.8)

    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(Adjusted p-value)')
    plt.title(f"Volcano Plot {cell_group_a} vs {cell_group_b}")
    plt.legend(handles=scatter.legend_elements()[0],
               labels=['Down', 'n.s.', 'Up'],
               title="Significance")
    if save_dir:
        if not os.path.exists(f"{save_dir}/{cell_group_a}_vs_{cell_group_b}"):
            os.makedirs(f"{save_dir}/{cell_group_a}_vs_{cell_group_b}")
        plt.savefig(f"{save_dir}/{cell_group_a}_vs_{cell_group_b}/volcano_plot.png", bbox_inches='tight',
                    dpi=300)
        print(f"{save_dir}/{cell_group_a}_vs_{cell_group_b}/volcano_plot.png successful saved")
        plt.close()
    else:
        plt.show()


def plot_merged_violin_comparison(
        adata: sc.AnnData,
        cell_group_a: str,
        cell_group_b: str,
        gene_lst: list,
        groupby: str = "cell_type",
        figsize: tuple = None,
        save_dir: str = None,
        dpi: int = 300
):
    """
    合并式分组比较图（多基因+双组别）

    参数：
        adata: AnnData对象
        cell_type_a: 细胞类型A
        cell_type_b: 细胞类型B
        gene_lst: 基因列表
        cell_type_key: 细胞类型列名
        figsize: 图片尺寸
        save_dir: 图片保存路径
        dpi: 分辨率

    返回：
        无返回，plot开头的函数只有绘图功能
    """
    import pandas as pd
    from statannotations.Annotator import Annotator
    import matplotlib.pyplot as plt
    import seaborn as sns
    # 数据校验
    valid_genes = [g for g in gene_lst if g in adata.var_names]
    if not valid_genes:
        raise ValueError("No valid genes found in adata.var_names")

    if save_dir is not None:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

    if figsize is None:
        figsize = tuple([0.8*len(valid_genes), 0.4*len(valid_genes)])

    # 构建绘图数据
    plot_data = []
    for gene in valid_genes:
        # 提取两种细胞类型的表达数据
        mask_a = (adata.obs[groupby] == cell_group_a)
        expr_a = adata[mask_a, gene].X.flatten()

        mask_b = (adata.obs[groupby] == cell_group_b)
        expr_b = adata[mask_b, gene].X.flatten()

        # 组装数据
        plot_data.extend([
            *[{'Gene': gene, 'Expression': e, 'Cell Type': cell_group_a} for e in expr_a],
            *[{'Gene': gene, 'Expression': e, 'Cell Type': cell_group_b} for e in expr_b]
        ])

    df = pd.DataFrame(plot_data)

    # 创建画布
    plt.figure(figsize=figsize)
    ax = plt.gca()

    # 绘制分组小提琴图
    sns.violinplot(
        x='Gene',
        y='Expression',
        hue='Cell Type',
        data=df,
        palette={cell_group_a: '#377eb8', cell_group_b: '#ff7f00'},
        split=True,  # 左右并列显示
        inner='quartile',  # 显示四分位线
        linewidth=1,
        ax=ax
    )

    # 关键修正部分：统计标注
    pairs = [((gene, cell_group_a), (gene, cell_group_b)) for gene in valid_genes]

    annotator = Annotator(
        ax=ax,
        pairs=pairs,
        data=df,
        x="Gene",
        y="Expression",
        hue="Cell Type",
        order=valid_genes  # 必须指定顺序
    )

    # 配置检验方法
    annotator.configure(
        test='Mann-Whitney',
        text_format='star',
        comparisons_correction="bonferroni"
    )

    # 应用标注
    annotator.apply_and_annotate()

    # 美化图形
    plt.title(f"Expression Comparison: {cell_group_a} vs {cell_group_b}")
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('')
    sns.despine()

    # 调整图例
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # 保存或显示
    if save_dir:
        if not os.path.exists(f"{save_dir}/{cell_group_a}_vs_{cell_group_b}"):
            os.makedirs(f"{save_dir}/{cell_group_a}_vs_{cell_group_b}")
        plt.savefig(f"{save_dir}/{cell_group_a}_vs_{cell_group_b}/merged_violin_comparison.png", bbox_inches='tight', dpi=dpi)
        print(f"{save_dir}/{cell_group_a}_vs_{cell_group_b}/merged_violin_comparison.png successful saved")
        plt.close()
    else:
        plt.show()


def cell_group_discrimination(cell_group_a="Primitive Streak", cell_group_b="Nascent mesoderm", groupby="cell_type", min_logfc=5.0,
        max_padj=0.01, save_dir="discrimination_plot"):
    """

    Do cell group gene expression comparison, the comparison is directed, with cell_group_a be the reference
    while cell_type_b be the target.
    The main method used in this function is scanpy.tl.rank_gene_groups.

    :param cell_group_a: cell group a grouped by the criteria groupby
    :param cell_group_b: cell group b grouped by the criteria groupby
    :param groupby: The criteria of grouping cells, such as stage, cell type.
    :param min_logfc: The minium logfoldchange, below this value will be filtered.
    :param max_padj: The maximum, bigger than this value will be filtered.
    :param save_dir: The dictionary where the figs to be stored, ALL csv files in this repo go directly to ./csv
    :return: None
    """
    # 寻找从Progenitor到Differentiated的差异基因
    transition_genes_df = find_transition_genes(
        adata,
        cell_group_a=cell_group_a,
        cell_group_b=cell_group_b,
        min_logfc=min_logfc,
        max_padj=max_padj,
        groupby=groupby,
    )

    # 输出结果
    print(f"Found {len(transition_genes_df)} transition genes")

    # plot_genes_df = pd.concat(objs=[transition_genes_df.sort_values("logfoldchanges", ascending=False).head(10),
    #                                 transition_genes_df.sort_values("pvals_adj", ascending=True).head(10)], axis=0)
    plot_genes_df = transition_genes_df.sort_values("pvals_adj", ascending=True).head(15)

    gene_lst = [
        "RTN3",
        "SMG1", "PANK3", "MCM7",
        "POLQ",
        "GALNT12", "RTN3", "TMEM87A", "COPB2",
        "FEN1",
        "MED1", "NCOA2",
        "SMG1", "CSNK2A1", "AKT2", "MKNK2", "TAOK2", "NEK7", "MAPKAPK2", "SIK3", "SIK2", "TAF1",
        "NDUFA7", "COX7B",
        "ATG3", "TOLLIP"]
    gene_lst = list(set(gene.capitalize() for gene in gene_lst))

    plot_merged_violin_comparison(adata=adata, gene_lst=plot_genes_df["names"].to_list(),
                                  groupby=groupby, cell_group_a=cell_group_a,
                                  cell_group_b=cell_group_b, save_dir=save_dir)

    plot_volcano(cell_group_a=cell_group_a, cell_group_b=cell_group_b, result_df=transition_genes_df, save_dir=save_dir)


def fate_decision_analysis(adata, parent_type, child1_type, child2_type,
                           n_genes=100, do_child_compare=False,**kwargs):
    """
    Enhanced version combining pseudotime analysis and visualization.

    1. Identifies genes correlated with fate decision using pseudotime
    2. Performs branch-specific differential expression
    3. Visualizes top candidate genes

    Pseudotime and PAGA trajectory will be recomputed in this function.
    """

    for i in [parent_type, child1_type, child2_type]:
        if i not in adata.obs['cell_type'].unique():
            raise ValueError(f"{i} not found in adata.obs['cell_type']. Please check params.")
        
    # Only keep three types of cells: [parent_type, child1_type, child2_type]
    sub_adata = adata[adata.obs["cell_type"].isin([parent_type, child1_type, child2_type]), :]

    def extract_diff_genes(f_adata, target, reference, f_n_genes):
        # 执行差异分析1
        sc.tl.rank_genes_groups(
            f_adata,
            groupby='cell_type',
            groups=[target],
            reference=reference,  # 直接比较两个子类型
            method='wilcoxon'
        )
        de_genes = sc.get.rank_genes_groups_df(f_adata, group=target)
        branch_genes = de_genes[de_genes["pvals"] < 0.001]['names'].tolist()
        # branch_genes = de_genes.head(f_n_genes)['names'].tolist()
        return branch_genes


    def save_txt(file_name, gene_lst):
        with open(file_name, 'w') as file:
            for item in gene_lst:
                file.write(f"{item}\n")


    # Find the genes that has small p-vals in both parent-child comparison.
    branch_genes_1 = extract_diff_genes(f_adata=sub_adata, target=child1_type, reference=parent_type, f_n_genes=n_genes)
    save_txt(f"{child1_type}_vs_{parent_type}.txt", branch_genes_1)
    branch_genes_2 = extract_diff_genes(f_adata=sub_adata, target=child2_type, reference=parent_type, f_n_genes=n_genes)
    save_txt(f"{child2_type}_vs_{parent_type}.txt", branch_genes_2)
    if do_child_compare:
        branch_genes_3 = extract_diff_genes(f_adata=sub_adata, target=child2_type, reference=child1_type, f_n_genes=n_genes)
        branch_genes_4 = extract_diff_genes(f_adata=sub_adata, target=child1_type, reference=child2_type, f_n_genes=n_genes)
        child_comp_genes = set(branch_genes_3) | set(branch_genes_4)
        candidate_genes = list(set(branch_genes_1) & set(branch_genes_2) & child_comp_genes)
        save_txt(f"{parent_type}_{child1_type}_{child2_type}.txt", candidate_genes)
    else:
        candidate_genes = list(set(branch_genes_1) & set(branch_genes_2))
    print(f"Found {len(candidate_genes)} significant genes in both two downstream directions.")


    # Visualize using matrixplot for clearer temporal patterns
    # sc.pl.matrixplot(
    #     sub_adata,
    #     var_names=candidate_genes,
    #     groupby='cell_type',
    #     standard_scale='var',
    #     cmap='viridis',
    #     title='Fate Decision Candidate Genes',
    #     save=".png"
    # )


if __name__ == "__main__":
    """
    本脚本用于比较在单细胞测序数据集中，对两类细胞群体之间的表达的基因进行显著性检验及绘图，找到两个群体之间的基因表达的差异。
    本脚本需要输入一个Anndata对象，其中包含单细胞测序数据，并且指定细胞按照何种方式分类，在该分类下的两个细胞群体的名称。
    The entry of this script is this module.
    """
    adata = sc.read("./h5ads/final_627_combat_re.h5ad")  # 加载预处理数据
    print(adata.n_vars)
    sc.settings.set_figure_params(dpi=300, dpi_save=300, facecolor="white")
    sc.settings.autosave = False
    sc.settings.autoshow = True
    sc.pp.filter_genes(adata, min_cells=300)
    # Epiblast Primitive Streak Nascent mesoderm
    # cell_group_discrimination(cell_group_a="Primitive Streak", cell_group_b="Anterior Primitive Streak",
    #                           max_padj=0.01, min_logfc=0.5,
    #                           groupby="cell_type", save_dir="discrimination_plot")
    fate_decision_analysis(adata, parent_type="Primitive Streak",
                           child1_type="Nascent mesoderm", child2_type="Anterior Primitive Streak", n_genes=400, do_child_compare=True)

