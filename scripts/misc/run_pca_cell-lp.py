import numpy as np
from scipy import sparse
from sklearn.decomposition import PCA
import umap
import pandas as pd
import os

data_file = 'all_sparse_matrices_10.txt'

data = pd.read_csv(data_file, delimiter = r"\s+", header=None)
mat_coo = sparse.coo_matrix((data[2], (data[0]-1, data[1])))
mat_coo = mat_coo.todense()

pca = PCA(n_components=1000, random_state=42)
pca.fit(mat_coo)
mat_transformed = pca.fit_transform(mat_coo)

out_mat_name = 'pca_coord_10.txt'
np.savetxt(out_mat_name, mat_transformed, fmt='%1.6f', delimiter="\t", newline="\n")

mat_transformed = mat_transformed[:,1:50]

reducer = umap.UMAP()
embedding = reducer.fit_transform(mat_transformed)
embedding.shape

out_umap_name = 'umap_coord_10_50.txt'
np.savetxt(out_umap_name, embedding, fmt='%1.6f', delimiter="\t", newline="\n")
