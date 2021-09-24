import sys
import scvi
import scanpy as sc
import h5py

# load data
adata = scvi.data.heart_cell_atlas_subsampled()

# print(adata.layers["counts"])
print(adata.X.copy())

# # check contents of data
f = h5py.File('data/hca_subsampled_20k.h5ad', 'r')
print('Keys: ')
print(list(f.keys()))

# # create dataset from 1 of the groups
fKey = f['X']
dset = fKey['data']

# # print ex dataset info
print("Dataset:")
print(dset)
print(dset.shape)
print(dset.dtype)
