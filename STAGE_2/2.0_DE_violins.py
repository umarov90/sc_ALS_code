import anndata as ad
import pandas as pd
from params import Params
import scanpy as sc
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
from gseapy.plot import gseaplot
from sanbomics.plots import volcano
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
matplotlib.use('macosx')


params = Params()
gene = "TAC3"
if params.deseq:
    res = pd.read_csv("/Users/ramzan/Dropbox/ALS_project/ALS_pycharm/STAGE_2/deseq_results.tsv",
                      sep="\t", index_col=0)
    sigs = res[(res.baseMean >= 10) & (res.padj < 0.05) ]
    # & (abs(res.log2FoldChange) > 0.5)
    dedf = sigs.sort_values('log2FoldChange', ascending = True).reset_index(drop = True)
else:
    df = pd.read_csv('/Users/ramzan/Dropbox/ALS_project/ALS_pycharm/R_scripts/diff_expr/output/all_cells'
                     '/glmQLFTest_result.tsv', sep="\t")
    significant_genes_edgeR = df[(df['logCPM'] > 5) &(df['FDR'] < 0.05) ]
    # & (abs(df['logFC']) > 0.5)
    dedf = significant_genes_edgeR.sort_values('logFC', ascending = True).reset_index(drop = True)

dedf = dedf[dedf['Symbol'] != 'MYC']

adata = ad.read_h5ad(params.file_path)
# adata = adata.raw.to_adata()
# subset = adata[adata.obs[params.cell_type] == "NSC"].copy()

temp = adata # adata[adata.obs[params.cell_type] == 'NSC']

i = np.where(temp.var_names == gene)[0][0]
a = temp[temp.obs.status == 'healthy'].X[:,i].todense()
b = temp[temp.obs.status == 'ALS'].X[:,i].todense()
p_value = stats.mannwhitneyu(np.asarray(a), np.asarray(b)).pvalue
fig, ax = plt.subplots(figsize=(7, 5))
sc.pl.violin(temp, gene, groupby='status', ax=ax, show=False)
ax.text(x=0.5, y=max(ax.get_ylim()),
        s=f'p-value: {p_value}', ha='center', va='bottom')
plt.savefig(gene + ".png", dpi=300)
# with open('ALS_GPT.txt') as f:
#     ALS_GPT = [x.strip() for x in list(f)]
#
# sc.tl.score_genes(subset, ALS_GPT, score_name = 'datp')
# sc.pl.violin(subset, 'datp', groupby='status')
# a = subset[subset.obs.status == 'healthy'].obs.datp.values
# b = subset[subset.obs.status == 'ALS'].obs.datp.values
# print(stats.mannwhitneyu(np.asarray(a), np.asarray(b)))
#
# sc.tl.pca(subset, use_highly_variable=False)
# sc.pp.neighbors(subset, n_neighbors=10, n_pcs=50)
# sc.tl.umap(subset)
# sc.pl.umap(subset, color = 'datp', vmax = 1)
