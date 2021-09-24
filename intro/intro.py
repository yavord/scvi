import sys
import scvi
import scanpy as sc


sc.set_figure_params(figsize=(4, 4))

# load data/check if its there
adata = scvi.data.heart_cell_atlas_subsampled()

print(adata)

# preprocessing: filter out anything with less than x counts
sc.pp.filter_genes(adata, min_counts=3)

adata.layers["counts"] = adata.X.copy() # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata # freeze the state in `.raw`

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="cell_source"
)
