from numpy import save
import scvi
import scanpy as sc

saveLocation = "models/intro/hca_ss"

# load data/check if its there
adata = scvi.data.heart_cell_atlas_subsampled()

print(adata.var.head())

# preprocessing: filter out anything with less than x counts
sc.pp.filter_genes(adata, min_counts=3)

# preserve current counts
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# save raw counts separately since normalized data is not used by scvi-tools models
adata.raw = adata

# select top 1200 genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="cell_source"
)

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    categorical_covariate_keys=["cell_source", "donor"],
    continuous_covariate_keys=["percent_mito", "percent_ribo"]
)

model = scvi.model.SCVI(adata)
# scvi.data.view_anndata_setup(model.adata)
model.train()
model.save(saveLocation, save_anndata=True)