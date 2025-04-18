# Final Year Project on Mouse Embryo Single Cell Analysis

------
README.md last updated on 18 Apr 2025. A lovely spring day.

The commit on 25 Apr 2025 deleted plenty outdated figures and folders.

---

## Repo Overview
This repo contains the code and main results of my final year proj. This proj is working on the single cell data of mouse embryo
from E6.5(i.e. day 6.5 after fertilization of mouse cells) to E7.5. The source data is downloaded from GEO with accession number GSE87038.

This proj has three goals:

1. Find mouse embryo cell trajectory.

2. Find the decisive genes which control the process. This part of work focus on the branching point of primitive streak cell.
Primitive streak cell has two downstream cell types, one is nascent mesoderm, the other is anterior primitive streak.
The former will develop into the mesoderm of the embryo, the latter will develop into the endoderm of the embryo.
Thus, primitive steak is a critical cell type in the trajectory.

3. Finally, after finding the decisive genes, do Gene Oncology analysis using NCBI David web tool and visualized the results

---

## File Description

### python scripts
A majority of function in this repo is annotated. 

***Many functions are not that versatile and may need to be carefully modified if they need to be reused in other proj.***

`sc_plot.py`: Python script contains plotting functions. **The function in this script are not used elsewhere.**

`pipeline.py`: **[The Main Entry Point of the Repo]** Python script contains single cell analysis pipeline, this script use many functions of 
`preprocess.py`, `umap_and_marker_gene.py`, `fa2_and_pseudotime.py`.

`preprocess.py`: Python script contains code on annotate cell, Gene ID Conversion, Gene ID transfer, sc data preprocess.
The function `preprocess` in `preprocess.py` applies the legacy [Scanpy Preprocess Workflow](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html).

`umap_and_marker_gene.py`: Python script contains code on clustering using umap, and find marker genes between different clusters
grouped by different criteria.

`fa2_and_pseudotime.py`: Python script contains code on clustering using PAGA and ForceAtlas2, and do pesudotime analysis
on ForceAtlas2 graph.

`compare_and_enrichment.py`: Python script contains code on finding and visualizing decisive genes.

`param_optim.py`: Python script contains code on adjusting scanpy.tl.umap hyper parameters to find the best fit ones.

### Result Folders

`{627|628}_{combat|mnn}_{re|no_re}`: Contains figures on different aspects.

{627|628}: Data range, 627 means using data from E6.5 to E7.5,
628 means using data from E6.5 to E8.5. 

{combat|mnn}: Batch correction method.

{re|no_re}: Do scanpy.regressout on sequencing depth(i.e. Anndata.obs['total_counts'])

`david_result`: Contains outputs from NCBI web tool David.

`insight_figures`: Contains output figures from `sc_plot.py`. Also, the default saving dir.

`discrimination plot`: Contains output figures from `compare_and_enrichment.py`. Also, the default saving dir.

### Pure Text Files

`README_2.md`: The second last version of `README.md`

`chat_with_deepseef.md`: Interesting insightful Q&A between me and DeepSeek V3 in Chinese.

----