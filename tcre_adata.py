import os
import pandas as pd
from scipy.io import mmread
import anndata as ad
from anndata import AnnData
from params import Params

params = Params()

samples_folder = "/osc-fs_home/ramzan/software/ALS/ALS_meta_tcre/"
adatas = []
for sample_folder in os.listdir(samples_folder):
    full_path = samples_folder + sample_folder + "/matrix/"
    if os.path.isdir(full_path):
        print(f'Processing files {sample_folder}')
        counts = mmread(full_path + 'matrix.mtx').T  # Transpose to make cells x genes
        genes = pd.read_csv(full_path + 'genes.tsv', header=None, sep='\t')[0].tolist()
        barcodes = pd.read_csv(full_path + 'barcodes.tsv', header=None, sep='\t')[0].tolist()
        obs = pd.read_csv(full_path + 'obs.csv', index_col=0)
        adata = AnnData(counts, obs=obs, var=pd.DataFrame(index=genes))
        adatas.append(adata)

print(f"Merging {len(adatas)} files...")
adata = ad.concat(adatas)
adata.obs['individual'] = adata.obs['individual'].astype(float).fillna(-1).astype(int).astype(str)
adata.write(params.folder + "tcre_adata.h5ad")
