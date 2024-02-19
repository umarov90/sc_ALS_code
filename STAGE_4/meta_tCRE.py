import os
import sys
import anndata as ad
import pandas as pd
from params import Params
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
import numpy as np
from scipy.io import mmread
import shutil

params = Params()
folder = "/Users/ramzan/als/tCRE_counts/"  # "/osc-fs_home/hon-chun/analysis/ramzan/ALS/scafe/count/merged_CRE_per_lib/"
out = "/Users/ramzan/als/tCRE_counts_out/"  # "/osc-fs_home/ramzan/software/ALS/ALS_meta_tcre/"
sample = "TFHS000444" # sys.argv[1]

adata = ad.read_h5ad(params.file_path)
adata = adata[adata.obs["sample_id"] == sample]

meta_all = adata.obs
meta_old_names = meta_all["old_obs_names"].tolist()
meta_new_names = adata.obs_names

matrix = mmread(folder + sample + "/matrix/matrix.mtx")
matrix = np.array(matrix.todense())
barcodes = pd.read_csv(folder + sample + "/matrix/barcodes.tsv", header=None, sep="\t")
barcodes.columns = ["barcodes"]
barcodes = barcodes["barcodes"].tolist()
barcodes = [s[:-2] for s in barcodes]
genes = pd.read_csv(folder + sample + "/matrix/genes.tsv", header=None, sep="\t")
genes.columns = ["start", "end"]
genes = genes["start"].tolist()
aggregated_matrices = []
for m in meta_old_names:
    m_list = m.split("\t")
    m_list = [s[len(sample) + 1:] for s in m_list]
    m_list = [barcodes.index(item) for item in m_list]
    aggregated_matrix = matrix[:, m_list].sum(axis=1)
    aggregated_matrices.append(aggregated_matrix)
a = np.vstack(aggregated_matrices).T
out_folder = f"{out}{sample}/matrix/"
os.makedirs(out_folder, exist_ok=True)
mmwrite(f'{out_folder}matrix.mtx', csr_matrix(a))
with open(f'{out_folder}barcodes.tsv', 'w') as file:
    for name in meta_new_names:
        file.write(name + '\n')
with open(f'{out_folder}barcodes_old.tsv', 'w') as file:
    for name in meta_old_names:
        file.write(name + '\n')
meta_all.to_csv(f'{out_folder}obs.csv')
shutil.copyfile(folder + sample + "/matrix/genes.tsv", out_folder + "genes.tsv")

