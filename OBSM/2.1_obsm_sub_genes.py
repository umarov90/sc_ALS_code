import anndata as ad
import numpy as np
import scanpy as sc
import scanpy.external as sce
from params import Params
import pandas as pd
import matplotlib
matplotlib.use('agg')

params = Params()
rn = 0
adata = ad.read_h5ad(params.file_path)
adata = adata[adata.obs["gender"] == "Female"].copy()
# print(len(adata))
# adata = adata[adata.obs["treatment"] != "Ropi"].copy()
# print(len(adata))
adata_small = adata[adata.obs["manual_anno_L0"].isin(['NSC', "pMN"])].copy() # | adata.obs["manual_anno_L0"].isin(['pMN'])
# adata_small = adata
# genes = pd.read_csv(p.folder + "classification_genes2.txt")["gene"].tolist()

if params.deseq:
    res = pd.read_csv("deseq_results.tsv", sep="\t", index_col=0)
else:
    df = pd.read_csv('/Users/ramzan/Dropbox/ALS_project/ALS_pycharm/R_scripts/diff_expr/output/all_cells'
                     '/glmQLFTest_result.tsv', sep="\t")
    res = df.rename(columns={"F":"stat"})
    res.index = res["Symbol"]

sigs = res[res.padj < 0.05].copy()
sigs['log2FoldChange'] = np.abs(sigs['log2FoldChange'])
sigs.sort_values('log2FoldChange', ascending=False, inplace=True)
genes = sigs.index.tolist()[:100]
genes.remove("MYC")
# genes = pd.read_csv(p.folder + "deseq_genes.tsv", sep="\t")["gene"].tolist()
# genes = list(set(genes + genes2))
adata_small = adata_small[:, genes]
# sc.tl.pca(adata, use_highly_variable=False, random_state=rn, n_comps=50)
sc.pp.neighbors(adata_small, n_neighbors=10, use_rep="X")
sc.tl.umap(adata_small, random_state=rn)
# sc.pl.umap(adata_small, color="day")
# sc.pl.umap(adata_small, color="pool")
sc.pl.umap(adata_small, color="manual_anno_L0", save="_NSC_pMN.pdf")
# sc.pl.umap(adata_small, color="manual_anno_L1")
sc.pl.umap(adata_small, color="status", save="_NSC_pMN_status.pdf")
# sc.pl.umap(adata_small, color="MYC", use_raw=False, save="_NSC_pMN_MYC.pdf")
# adata.obsm['X_umap_deseq'] = adata_small.obsm['X_umap']
# adata.obsm['X_umap_status_classification'] = np.zeros((len(adata), 2), dtype=np.float64)
# adata.obsm['X_umap_status_classification'][adata.obs["manual_anno_L0"].isin(['iPSC', 'pMN', 'NSC', 'oMN'])] = adata_small.obsm['X_umap']
# cm.safe_save(adata, p.file_path)
