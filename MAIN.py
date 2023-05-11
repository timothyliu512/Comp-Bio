import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

results_file = 'write/Analysis_Results'  # the file that will store the analysis results

adata = sc.read_10x_mtx(
    'data',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

adata.var_names_make_unique() 

# Preprocessing
adata.raw = adata
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['n_genes'] < 2500, :]
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var['highly_variable']]

# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# Louvain clustering
sc.tl.louvain(adata)
sc.tl.umap(adata)

sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Visualizing the clustering
sc.pl.umap(adata, color='louvain', title='UMAP visualization of Louvain clusters', legend_loc='on data')

adata

# Extract top ranked genes for each cluster
marker_genes = {}
for cluster in range(len(adata.obs['louvain'].cat.categories)):
    genes_of_interest = adata.uns['rank_genes_groups']['names'][str(cluster)]
    pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][str(cluster)]
    
    # Filter genes with adjusted p-value < 0.05
    genes_significant = genes_of_interest[pvals_adj < 0.05]
    
    # Add to dictionary
    marker_genes[cluster] = genes_significant
    
# Print top 10 marker genes for each cluster
for cluster, genes in marker_genes.items():
    print(f"Cluster {cluster}:")
    print(genes[:10])
    print("\n")

