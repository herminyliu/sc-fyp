import os

import scanpy as sc
import pandas as pd

import trajectory_inference


def find_transition_genes(
        adata: sc.AnnData,
        cell_type_a: str,
        cell_type_b: str,
        groupby: str = "cell_type",
        min_logfc: float = 0.25,
        max_padj: float = 0.05,
        key_added: str = "transition_genes"
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

    返回：
        包含显著差异基因及统计量的DataFrame
    """
    # 校验输入
    if groupby not in adata.obs:
        raise ValueError(f"Column '{groupby}' not found in adata.obs_names:{adata.obs}")

    cell_types = adata.obs[groupby].unique()
    if cell_type_a not in cell_types:
        raise ValueError(f"Cell type '{cell_type_a}' not found")
    if cell_type_b not in cell_types:
        raise ValueError(f"Cell type '{cell_type_b}' not found")

    # 创建细胞子集
    mask = adata.obs[groupby].isin([cell_type_a, cell_type_b])
    adata_sub = adata[mask].copy()

    # 设置比较组
    adata_sub.obs["comparison_group"] = [
        "target" if ct == cell_type_b else "reference"
        for ct in adata_sub.obs[groupby]
    ]

    # 进行差异表达分析
    sc.tl.rank_genes_groups(
        adata_sub,
        groupby="comparison_group",
        groups=["target"],
        reference="reference",
        method="wilcoxon",
        use_raw=False  # 使用已经归一化后的数据
    )

    # 提取结果
    result = sc.get.rank_genes_groups_df(adata_sub, group="target")

    # 过滤显著基因
    # significant_genes = result[
    #     (result["logfoldchanges"] >= min_logfc) &
    #     (result["pvals_adj"] <= max_padj)
    #     ].sort_values("logfoldchanges", ascending=False)
    significant_genes = result.sort_values(by="logfoldchanges", ascending=False)

    # 存储结果到原始adata对象
    adata.uns[key_added] = significant_genes

    return significant_genes


def plot_volcano(result_df,
                 logfc_thresh=0.5,
                 padj_thresh=0.05,
                 title="Volcano Plot"):
    """
    绘制差异基因火山图
    """
    import matplotlib.pyplot as plt
    from numpy import log10
    plt.figure(figsize=(8, 6))

    # 创建颜色分类列
    result_df['significance'] = 'n.s.'
    result_df.loc[(result_df['logfoldchanges'] >= logfc_thresh) &
                  (result_df['pvals_adj'] <= padj_thresh), 'significance'] = 'Up'
    result_df.loc[(result_df['logfoldchanges'] <= -logfc_thresh) &
                  (result_df['pvals_adj'] <= padj_thresh), 'significance'] = 'Down'

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
    plt.title(title)
    plt.legend(handles=scatter.legend_elements()[0],
               labels=['Down', 'n.s.', 'Up'],
               title="Significance")
    plt.show()


def plot_merged_comparison(
        adata: sc.AnnData,
        cell_type_a: str,
        cell_type_b: str,
        gene_lst: list,
        cell_type_key: str = "cell_type",
        figsize: tuple = None,
        save_path: str = None,
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
        save_path: 保存路径
        dpi: 分辨率
    """
    import pandas as pd
    from statannotations.Annotator import Annotator
    import matplotlib.pyplot as plt
    import seaborn as sns
    # 数据校验
    valid_genes = [g for g in gene_lst if g in adata.var_names]
    if not valid_genes:
        raise ValueError("No valid genes found in adata.var_names")

    if save_path is not None:
        if not os.path.exists(save_path):
            os.makedirs(save_path)

    if figsize is None:
        figsize = tuple([1.2*len(valid_genes), 0.6*len(valid_genes)])

    # 构建绘图数据
    plot_data = []
    for gene in valid_genes:
        # 提取两种细胞类型的表达数据
        mask_a = (adata.obs[cell_type_key] == cell_type_a)
        expr_a = adata[mask_a, gene].X.flatten()

        mask_b = (adata.obs[cell_type_key] == cell_type_b)
        expr_b = adata[mask_b, gene].X.flatten()

        # 组装数据
        plot_data.extend([
            *[{'Gene': gene, 'Expression': e, 'Cell Type': cell_type_a} for e in expr_a],
            *[{'Gene': gene, 'Expression': e, 'Cell Type': cell_type_b} for e in expr_b]
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
        palette={cell_type_a: '#377eb8', cell_type_b: '#ff7f00'},
        split=True,  # 左右并列显示
        inner='quartile',  # 显示四分位线
        linewidth=1,
        ax=ax
    )

    # 关键修正部分：统计标注
    pairs = [((gene, cell_type_a), (gene, cell_type_b)) for gene in valid_genes]

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
    plt.title(f"Expression Comparison: {cell_type_a} vs {cell_type_b}")
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('')
    sns.despine()

    # 调整图例
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # 保存或显示
    if save_path:
        plt.savefig(f"{save_path}/merged_comparison.png", bbox_inches='tight', dpi=dpi)
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    # 示例调用
    adata = sc.read("final_627_combat_re.h5ad")  # 加载预处理数据
    sc.settings.set_figure_params(dpi=300, dpi_save=300, facecolor="white")
    sc.pp.filter_genes(adata, min_cells=300)

    # 寻找从Progenitor到Differentiated的差异基因
    transition_genes = find_transition_genes(
        adata,
        cell_type_a="Primitive Streak",
        cell_type_b="Nascent mesoderm",
        min_logfc=0.5,
        max_padj=0.01,
        groupby="cell_type",
    )

    # 输出结果
    print(f"Found {len(transition_genes)} transition genes")
    significant_genes = transition_genes[
        (transition_genes["logfoldchanges"] >= 0.25) &
        (transition_genes["pvals_adj"] <= 0.05)]
    print(f"Found {len(significant_genes)} significant transition genes")

    # plot_volcano(result_df=transition_genes)

    plot_genes_df = pd.concat(objs=[transition_genes.sort_values("logfoldchanges", ascending=False).head(10),
                 transition_genes.sort_values("pvals_adj", ascending=True).head(10)], axis=0)

    # plot_genes_df = transition_genes.sort_values("logfoldchanges", ascending=False).head(10)
    # trajectory_inference.plot_pseudotime(adata, plot_genes_df["names"].to_list(), cell_sample_step=1)
    # plot_genes_df = transition_genes.sort_values("pvals_adj", ascending=True).head(15)
    # trajectory_inference.plot_pseudotime(adata, plot_genes_df["names"].to_list(), cell_sample_step=1)
    plot_merged_comparison(adata=adata, gene_lst=plot_genes_df["names"].to_list(),
                                    cell_type_key="cell_type", cell_type_a="Primitive Streak",
                                    cell_type_b="Nascent mesoderm", save_path="voline_plots")
