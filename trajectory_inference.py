import scanpy as sc

if __name__ == "__main__":
    adata = sc.read_h5ad('clustered_data_123_noBatchEffect.h5ad')
    sc.set_figure_params(dpi=150)

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
    groupby = "cell_time"
    # paga图，图元素包含有抽象的点和边
    sc.tl.paga(adata, groups=groupby)
    sc.pl.paga(adata, color=[groupby, "Camk1g", "Abca1", "Gm16120"])

    # paga散点图，图元素包含很多散点
    sc.tl.draw_graph(adata, init_pos="paga")
    sc.pl.draw_graph(
        adata, color=[groupby, "Camk1g", "Abca1", "Gm16120"], legend_loc="on data"
    )
