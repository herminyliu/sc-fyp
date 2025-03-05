import scanpy as sc

adata = sc.read_h5ad('clustered_data_123_noBatchEffect.h5ad')
print(adata.obs.n_genes[1])
print(adata.obs.n_genes_by_counts[1])
print(adata.obs.total_counts[1])
print(adata)