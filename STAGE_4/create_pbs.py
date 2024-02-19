import anndata as ad
import pandas as pd
from params import Params
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.plots import volcano
import numpy as np
import matplotlib

params = Params()

adata = ad.read_h5ad(params.folder + "tcre_adata.h5ad")
column_name_level = "manual_anno_L2"
for cluster in adata.obs[column_name_level].unique():
    print(cluster)
print(len(adata.obs[column_name_level].unique()))
adata.obs['individual'] = adata.obs['individual'].astype(float).fillna(-1).astype(int).astype(str)

pbs = []
for cluster in adata.obs[column_name_level].unique():
    samp_cell_subset = adata[adata.obs[column_name_level] == cluster]

    rep_adata = sc.AnnData(X=samp_cell_subset.X.sum(axis=0),
                           var=samp_cell_subset.var[[]])
    rep_adata.obs_names = column_name_level
    pbs.append(rep_adata)
pb = sc.concat(pbs)

os.makedirs(out_folder, exist_ok=True)
mmwrite(f'{out_folder}matrix.mtx', csr_matrix(a))

# genes
with open(f'{out_folder}barcodes.tsv', 'w') as file:
    for name in meta_new_names:
        file.write(name + '\n')