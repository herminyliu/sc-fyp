import numpy as np
import pandas as pd
from scipy.sparse import issparse
from sklearn.linear_model import LinearRegression
import scanpy as sc


def select_genes_by_dpt_trend(adata, n_genes=30):
    """
    Select genes based on their expression trend along dpt_order_indices.

    Parameters:
        adata (AnnData): Annotated data matrix.
        n_genes (int): Number of top genes to select (default: 10).

    Returns:
        selected_genes (list): List of selected gene names.
    """
    # Check if dpt_order_indices is available
    if "dpt_order_indices" not in adata.obs:
        raise ValueError("dpt_order_indices not found in adata.obs. Please compute it first!")

    # Check if highly variable genes are already computed
    if "highly_variable" not in adata.var:
        raise ValueError("Highly variable genes not found in adata.var. Please compute them first!")

    # Subset to highly variable genes
    adata_hvg = adata[:, adata.var.highly_variable]

    # Subset to genes contain in over 5000 cells, more strict requirement.
    sc.pp.filter_genes(adata_hvg, min_cells=20)

    # Get dpt_order_indices and expression matrix
    dpt_order = adata_hvg.obs["dpt_order_indices"].values
    X = adata_hvg.X.toarray() if issparse(adata_hvg.X) else adata_hvg.X  # Ensure dense matrix

    # Calculate the slope (trend) for each gene
    slopes = []
    for i in range(X.shape[1]):
        y = X[:, i]  # Expression values of the current gene
        model = LinearRegression()
        model.fit(dpt_order.reshape(-1, 1), y)  # Fit linear regression
        slopes.append(model.coef_[0]) # Slope

    # Map slopes to gene names
    gene_slopes = pd.DataFrame({"gene": adata_hvg.var_names, "slope": slopes})

    # Select top genes
    # Positive slope: increasing expression; Negative slope: decreasing expression
    top_genes_positive = gene_slopes.nlargest(n_genes, "slope")  # Top n_genes with increasing expression
    top_genes_negative = gene_slopes.nsmallest(n_genes, "slope")  # Top n_genes with decreasing expression

    # Combine results
    top_genes = pd.concat([top_genes_positive, top_genes_negative])

    # Print results
    print("Selected genes based on expression trend:")
    print(top_genes)

    # Extract gene names
    selected_genes = top_genes["gene"].tolist()
    return selected_genes
