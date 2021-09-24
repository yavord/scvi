import scvi
import scanpy as sc


pbmc3k = scvi.data.read_h5ad("data/pbmc3k_raw.h5ad")
pbmc5k = sc.read_10x_h5(
    "data/5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5",
    gex_only=False
)

print('pbmc3k:')
print(pbmc3k)
print(pbmc3k.var.head())

print('pbmc5k:')
print(pbmc5k)
print(pbmc5k.var.head())
