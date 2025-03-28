import preprocess
import umap_clustering
import scanpy as sc


def umap_param_optim():
    adata = sc.read_h5ad("raw_without_batch_combat_no_re.h5ad")

    sc.settings.figdir = "umap_optim_results"
    sc.settings.autoshow = False
    sc.settings.autosave = True
    for n_neighbors in range(10, 25, 5):
        for min_dist in range(1, 9, 2):
            umap_clustering.cluster_umap(u_adata=adata, leiden_resolution=0.9, n_neighbors=n_neighbors, spread=1.0, min_dist=min_dist/10)
            umap_clustering.plot_umap(pp_adata=adata, save_fig_path=f"{n_neighbors}_{min_dist/10}_umap.png")




if __name__ == "__main__":
    umap_param_optim()