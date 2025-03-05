import scanpy as sc
from anndata import AnnData

import preprocess


if __name__ == '__main__':
    sc.settings.figdir = "./show_batch_effect"
    sc.settings.autosave = True
    sc.settings.autoshow = False
    adata = sc.read_h5ad("clustered_data_123.h5ad")
    adata = preprocess.annotate_cells(adata)
    sc.pl.umap(adata=adata, color=["cell_type", "cell_time", "batch", "leiden"])
