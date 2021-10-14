import scvi
import scanpy as sc

saveLocation = "models/intro/"
sc.set_figure_params(figsize=(4, 4))

adata = scvi.data.read_h5ad(saveLocation+"adata.h5ad")
model = scvi.model.SCVI.load(saveLocation, adata = adata)

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

# print(adata[adata.obs.cell_type])
adata_subset = adata[adata.obs.cell_type == "Fibroblast"]
latent_subset = model.get_latent_representation(adata_subset)

denoised = model.get_normalized_expression(adata_subset, library_size=1e4)
# print(denoised.iloc[:5, :5])
adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=10e4)


### UMAP without batch correction
# run PCA -> UMAP plots
# sc.tl.pca(adata, svd_solver='arpack')
# sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
# sc.tl.umap(adata, min_dist=0.3)

# sc.pl.umap(
#     adata,
#     color=["cell_type"],
#     frameon=False,
# )
# sc.pl.umap(
#     adata,
#     color=["donor", "cell_source"],
#     ncols=2,
#     frameon=False,
# )

### UMAP with batch correction
sc.pp.neighbors(adata, use_rep="X_scVI") # use latent space
sc.tl.umap(adata, min_dist=0.3)

sc.pl.umap(
    adata,
    color=["cell_type"],
    frameon=False,
)
sc.pl.umap(
    adata,
    color=["donor", "cell_source"],
    ncols=2,
    frameon=False,
)