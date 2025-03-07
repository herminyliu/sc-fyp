import scanpy as sc

FILE_PATH = "clustered_data_123_noBatchEffect.h5ad"
SAVING_FIG_FOLDER = './preprocessed_figures_123_noBatchEffect'
IS_AUTOSAVE = False
IS_AUTOSHOW = True


def paga(p_adata, color, groups):
    # paga图，图元素包含有抽象的点和边
    # sc.tl.paga不包含groupby参数
    sc.tl.paga(p_adata, groups=groups)
    # sc.pl.paga不包含groups和groupby这两种参数。groups参数自动继承自上面的tl.paga，自动着色出分组。
    sc.pl.paga(
        p_adata, color=color
    )
    return p_adata


def paga_scatter(ps_adata, color):
    # paga散点图，图元素包含很多散点
    # sc.tl.draw_graph不包含groups和groupby这两种参数。groups参数自动继承自前面的tl.paga
    sc.tl.draw_graph(ps_adata, init_pos="paga")
    # sc.pl.draw_graph不包含groups和groupby这两种参数
    sc.pl.draw_graph(
        adata=ps_adata, color=color, legend_loc="on data"
    )
    return ps_adata


def setting():
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
    # setting()
    adata = sc.read_h5ad(FILE_PATH)
    print()

    # 对细胞类型进行缩写处理，避免图上元素覆盖交叠
    adata.obs["cell_type_abbr"] = adata.obs["cell_type"].cat.rename_categories(
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
    # 对leiden标签进行注释，0-5可以归为一类
    adata.obs["leiden_anno"] = adata.obs["leiden"].cat.rename_categories(
        {
            "0": "0/PS",
            "1": "1/PS",
            "2": "2/PS",
            "3": "3/PS",
            "4": "4/PS",
            "5": "5/PS",
        }
    )
    groups = "cell_type_abbr"
    marker_genes = [groups, "Fgf3", "Gm17416", "Gm43604"]
    # 是否一定要先运行paga才再能运行paga_scatter呢？sc.tl.paga, sc.pl.paga, sc.tl.draw_graph, sc.pl.draw_graph函数逻辑复杂。
    # 实验证明不需要！sc.tl.draw_graph并不依赖sc.tl.paga。sc.pl.draw_graph会自动判断传入的color参数是否为var中的基因，或对obs的注释。
    # TODO：
    # 如果为obs的categories类注释，则会自动细胞群体着色。如果为var_name中的基因，那么自动按照各单细胞中该基因的表达量着色。那么着色的值是什么？是scale后的那个表达量吗？
    # 而且draw_graph的得到的散点图的横轴纵轴的含义是？散点的分布和UMAP图的分布明显不一样，那么这又是什么算法？
    # adata = paga(p_adata = adata, color=marker_genes, groups=groups)
    # 先过周末吧
    adata = paga_scatter(ps_adata=adata, color=marker_genes)
