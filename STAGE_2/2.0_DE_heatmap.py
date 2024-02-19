import anndata as ad
import pandas as pd
from params import Params
import scanpy as sc
import numpy as np
import matplotlib
matplotlib.use('macosx')


params = Params()


if params.deseq:
    res = pd.read_csv("/Users/ramzan/Dropbox/ALS_project/ALS_pycharm/STAGE_2/deseq_results.tsv",
                      sep="\t", index_col=0)
    sigs = res[(res.padj < 0.05)]
    # (res.baseMean >= 1) & (abs(res.log2FoldChange) > 0.5) &
    dedf = sigs.sort_values('log2FoldChange', ascending=False).reset_index(drop=True)
else:
    df = pd.read_csv('/Users/ramzan/Dropbox/ALS_project/ALS_pycharm/R_scripts/diff_expr/output/all_cells'
                     '/glmQLFTest_result.tsv', sep="\t")
    sigs = df[(df['logCPM'] > 5) & (df['FDR'] < 0.05) ]
    # & (abs(df['logFC']) > 0.5)
    dedf = sigs.sort_values('logFC', ascending=True).reset_index(drop=True)
print(len(sigs))
# dedf = dedf[dedf['Symbol'] != 'MYC']
adata = ad.read_h5ad(params.file_path)
adata = adata.raw.to_adata()
adata = adata[adata.obs["gender"] == "Female"].copy()
adata = adata[adata.obs[params.cell_type].isin(["NSC", "pMN"])].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

groupby = ["treatment", "status"][1]
np.unique(adata.obs[groupby])
if groupby == "treatment":
    vals = ["Ctrl", "Ropi"]
    print(len(adata))
    adata = adata[adata.obs["status"] == "ALS"].copy()
    adata = adata[adata.obs["day"] != "D0"].copy()
    print(len(adata))
elif groupby == "status":
    vals = ["healthy", "ALS"]
    print(len(adata))
    adata = adata[adata.obs["treatment"] != "Ropi"].copy()
    print(len(adata))

def print_expression(gene_name):
    i = np.where(adata.var_names == gene_name)[0][0]
    a = adata[adata.obs[groupby] == vals[0]].X[:, i]
    b = adata[adata.obs[groupby] == vals[1]].X[:, i]
    print(f"{gene_name} expression: {vals[0]}: {a.mean()} {vals[1]}: {b.mean()}")


most_up = dedf.iloc[0].Symbol
print(print_expression(most_up))
top_n = 20
genes_to_show = dedf[-top_n:].Symbol.tolist() + dedf[:top_n].Symbol.tolist()
# for g in genes_to_show:
#     print_expression(g)

sc.pl.heatmap(adata, genes_to_show, groupby=groupby,
              swap_axes=True, save=f"top{top_n}_de_genes.pdf", use_raw=False)

adata = ad.read_h5ad(params.file_path)
adata = adata.raw.to_adata()
adata = adata[adata.obs["gender"] == "Female"].copy()
adata = adata[adata.obs[params.cell_type].isin(["NSC", "pMN"])].copy()
pbs = []
for i in adata.obs["individual"].unique():
    samp_cell_subset = adata[adata.obs["individual"] == i]
    rep_adata = sc.AnnData(X=samp_cell_subset.X.sum(axis=0),
                           var=samp_cell_subset.var[[]])
    rep_adata.obs_names = [i]
    rep_adata.obs[groupby] = samp_cell_subset.obs[groupby].iloc[0]
    pbs.append(rep_adata)
pb = sc.concat(pbs)
pb.X = np.array(pb.X)
sc.pp.normalize_total(pb, target_sum=1e4)
sc.pp.log1p(pb)
sc.pp.scale(pb, max_value=10)
sc.pl.heatmap(pb, genes_to_show, groupby=groupby,
              swap_axes=True, save=f"top{top_n}_de_genes_pb.pdf", use_raw=False)


pbs = []
for g in adata.obs[groupby].unique():
    samp_cell_subset = adata[adata.obs[groupby] == g]
    rep_adata = sc.AnnData(X=samp_cell_subset.X.sum(axis=0),
                           var=samp_cell_subset.var[[]])
    rep_adata.obs_names = [g]
    rep_adata.obs[groupby] = samp_cell_subset.obs[groupby].iloc[0]
    pbs.append(rep_adata)
pb = sc.concat(pbs)
pb.X = np.array(pb.X)
# sc.pp.normalize_total(pb, target_sum=1e4)
sc.pp.log1p(pb)
sc.pl.heatmap(pb, genes_to_show, groupby=groupby,
              swap_axes=True, save=f"top{top_n}_de_genes_pb.pdf", use_raw=False)