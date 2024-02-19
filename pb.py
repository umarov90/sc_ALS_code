import anndata as ad
import pandas as pd
from params import Params
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix

params = Params()

adata = ad.read_h5ad(params.folder + "tcre_adata.h5ad")
column_name_level = "manual_anno_L2"
for cluster in adata.obs[column_name_level].unique():
    print(cluster)
print(len(adata.obs[column_name_level].unique()))
exit()
pbs = []
for individual in adata.obs["individual"].unique():
    adata_i = adata[adata.obs["individual"] == individual]
    for cluster in adata_i.obs[column_name_level].unique():
        samp_cell_subset = adata_i[adata_i.obs[column_name_level] == cluster]
        rep_adata = sc.AnnData(X=samp_cell_subset.X.sum(axis=0),
                               var=samp_cell_subset.var[[]])
        rep_adata.obs_names = [individual + "_" + cluster]
        rep_adata.obs['individual'] = samp_cell_subset.obs['individual'].iloc[0]
        rep_adata.obs[column_name_level] = samp_cell_subset.obs[column_name_level].iloc[0]
        pbs.append(rep_adata)
print(len(pbs))
pb = sc.concat(pbs)
pb.X = csr_matrix(pb.X)
print(np.max(pb.X))
print(pb.shape)
print(pb.X)
print(pb.X.dtype)
pb.X = pb.X.astype('float64')
pb.write(params.folder + "individual_manual_anno_L2_counts.h5ad")
sc.pp.normalize_per_cell(pb, counts_per_cell_after=1e6)
pb.write(params.folder + "individual_manual_anno_L2_norm_myself.h5ad")