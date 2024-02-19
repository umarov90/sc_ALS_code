import anndata as ad
import pandas as pd
from params import Params
import gseapy as gp
from gseapy.plot import gseaplot
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')
from matplotlib import cm
import seaborn as sns
import textwrap
import numpy as np

params = Params()

if params.deseq:
    res = pd.read_csv("deseq_results.tsv",
                      sep="\t", index_col=0)
else:
    df = pd.read_csv('/Users/ramzan/Dropbox/ALS_project/ALS_pycharm/R_scripts/diff_expr/output/all_cells'
                     '/glmQLFTest_result.tsv', sep="\t")
    res = df.rename(columns={"F":"stat"})
    res.index = res["Symbol"]
res = res[(res.padj < 0.05)]

# GSEA
ranking = res[['Symbol', 'log2FoldChange']].dropna()
ranking['log2FoldChange'] = -1 * np.abs(ranking['log2FoldChange'])
ranking = ranking.drop_duplicates('Symbol')
ranking.index = ranking.index.str.upper()
pre_res = gp.prerank(rnk=ranking.drop(columns=['Symbol']),
                         gene_sets="../data/GMTs/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
# print(pre_res.res2d[pre_res.res2d["Term"] == "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS"]["NES"].iloc[0])

out = []

for term in list(pre_res.results):
    out.append([term,
                pre_res.results[term]['fdr'],
                pre_res.results[term]['es'],
                pre_res.results[term]['nes']])


out_df = pd.DataFrame(out, columns=['term', 'fdr', 'es', 'nes']).sort_values('nes', ascending=False).reset_index(drop=True)
term = "KEGG_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS"
print(out_df[out_df["term"]==term])
gseaplot(rank_metric=pre_res.ranking, figsize=(6, 10), ofname=f"GSEA_ALS_KEGG.png",
         term=term,
         **pre_res.results[term])
out_df = out_df[out_df['fdr'] < 0.5]
out_df = out_df.head(50)
fig, (ax, cax) = plt.subplots(figsize=(12, 18), ncols=2, gridspec_kw={'width_ratios': [20, 1], 'wspace': 0.05})
cmap = mpl.cm.bwr_r
norm = mpl.colors.Normalize(vmin=out_df.fdr.min(), vmax=out_df.fdr.max())
mapper = cm.ScalarMappable(norm=norm, cmap=cm.bwr_r)
sns.barplot(data=out_df, x='nes', y='term', palette=mapper.to_rgba(out_df.fdr.values), ax=ax)
ax.set_yticklabels([textwrap.fill(e, 30) for e in out_df['term']])
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
cb.set_label('FDR')
plt.subplots_adjust(left=0.4)
plt.tight_layout()
plt.savefig("GSEA.png", dpi=300)
