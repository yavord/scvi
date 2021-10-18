import scvi
import scanpy as sc

saveLocation = "models/intro/pbmc"

# load/view pbmc3k
pbmc3k = scvi.data.read_h5ad("data/pbmc3k_raw.h5ad")
# print('pbmc3k:')
# print(pbmc3k)
# print(pbmc3k.var.head())

# load/view pbmc5k
pbmc5k = sc.read_10x_h5(
    "data/5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5",
    gex_only=False
)
pbmc5k.var_names_make_unique()
scvi.data.organize_cite_seq_10x(pbmc5k)

# print('pbmc5k:')
# print(pbmc5k)
# print(pbmc5k.var.head())

# combine both datasets
adata = pbmc5k.concatenate(pbmc3k)
# print(adata.obs.sample(n=5)) # check batch keys

# filter out < 3 counts 
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_counts=3)

# keep current counts
adata.layers["counts"] = adata.X.copy()

# normalize data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# keep raw data
adata.raw = adata

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
print(scvi.data.view_anndata_setup(adata))
# scvi.model.TOTALVI.setup_anndata(pbmc5k, protein_expression_obsm_key="protein_expression")
model = scvi.model.SCVI(adata)
model.train()
model.save(saveLocation, save_anndata=True)

