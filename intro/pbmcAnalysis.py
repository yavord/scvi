import scvi
import scanpy as sc

saveLocation = "models/intro/pbmc/"
sc.set_figure_params(figsize=(4, 4))

adata = scvi.data.read_h5ad(saveLocation+"adata.h5ad")
model = scvi.model.SCVI.load(saveLocation, adata = adata)

latent = model.get_latent_representation()
adata.obsm["X_scVI"] = latent

denoised = model.get_normalized_expression(adata, library_size=1e4)
adata.layers["scvi_normalized"] = model.get_normalized_expression(library_size=10e4)

sc.pp.neighbors(adata, use_rep="X_scVI") # use latent space
sc.tl.leiden(adata)
sc.tl.umap(adata, min_dist=0.3)

sc.pl.umap(
    adata,
    color=["batch", "leiden"],
    frameon=False,
    ncols = 1,
)

de_df = model.differential_expression(
    groupby="batch"
)

markers = {}
cats = adata.obs.batch.cat.categories

for i in cats:
    cid = "{} vs Rest".format(i)
    batch_df = de_df.loc[de_df.comparison == cid]

    batch_df = batch_df[batch_df.lfc_mean > 0]

    batch_df = batch_df[batch_df["bayes_factor"] > 3]
    batch_df = batch_df[batch_df["non_zeros_proportion1"] > 0.1]

    markers[i] = batch_df.index.tolist()[:3]

sc.pl.dotplot(
    adata,
    markers,
    groupby='batch',
    dendrogram=True,
    color_map="Blues",
    swap_axes=True,
    use_raw=True,
    standard_scale="var",
)

sc.pl.heatmap(
    adata,
    markers,
    groupby='batch',
    layer="scvi_normalized",
    standard_scale="var",
    dendrogram=True,
    figsize=(8, 12)
)
