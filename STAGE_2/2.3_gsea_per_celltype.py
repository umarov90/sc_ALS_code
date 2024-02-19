import os
import sys
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import gseapy
import traceback
import joblib
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.plots import volcano
from gseapy.plot import gseaplot
from scipy import stats
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import textwrap
script_dir = os.path.dirname(__file__)  # Directory containing myscript.py
parent_dir = os.path.dirname(script_dir)  # Project directory
sys.path.append(parent_dir)
import utils.common as cmf
from utils.params import Params


params = Params()
adata = ad.read_h5ad(params.file_path)
adata = adata.raw.to_adata()
column_name_level = "manual_anno_L2"

groupby = ["treatment", "status"][0]
iter_obs = ["individual", params.cell_type][1]
fig_folder = f"per_{iter_obs}_{groupby}_DE/"
np.unique(adata.obs[groupby])
if groupby == "treatment":
    vals = ["Ctrl", "Ropi"]
    print(len(adata))
    adata = adata[adata.obs["status"] == "ALS"].copy()
    adata = adata[adata.obs["day"].isin(["D14"])].copy()
    print(len(adata))
elif groupby == "status":
    vals = ["healthy", "ALS"]
    print(len(adata))
    adata = adata[adata.obs["treatment"] != "Ropi"].copy()
    print(len(adata))

adata.obs['individual'] = adata.obs['individual'].astype(str) + '_' + \
                          adata.obs['gender'].astype(str).str[0] + '_' + \
                          adata.obs['status'].astype(str).apply(lambda x: 'H' if x == 'healthy' else 'A')

sig_size = {}
genes_heatmap = {}
enr_results = []
for cluster in list(adata.obs[iter_obs].unique()) + ["everything"]:
    print(cluster)
    if cluster == "everything":
        adata_c = adata
    else:
        adata_c = adata[adata.obs[iter_obs] == cluster].copy()
    try:
        pbs = []
        for individual in adata_c.obs["individual"].unique():
            adata_i = adata_c[adata_c.obs['individual'] == individual]
            for subcluster in adata_i.obs[column_name_level].unique():
                adata_s = adata_i[adata_i.obs[column_name_level] == subcluster]
                for gr in adata_s.obs[groupby].unique():
                    samp_cell_subset = adata_s[adata_s.obs[groupby] == gr]
                    rep_adata = sc.AnnData(X=samp_cell_subset.X.sum(axis=0),
                                           var=samp_cell_subset.var[[]])
                    rep_adata.obs_names = [individual + "_" + subcluster + "_" + gr]
                    rep_adata.obs[groupby] = samp_cell_subset.obs[groupby].iloc[0]
                    rep_adata.obs['gender'] = samp_cell_subset.obs['gender'].iloc[0]
                    rep_adata.obs['manual_anno_L2'] = samp_cell_subset.obs['manual_anno_L2'].iloc[0]
                    pbs.append(rep_adata)
        pb = sc.concat(pbs)
        counts = pd.DataFrame(pb.X, columns=pb.var_names)
        design_factors = [groupby, "manual_anno_L2"]
        if iter_obs == params.cell_type:
            design_factors.append("gender")
        dds = DeseqDataSet(
            counts=counts,
            metadata=pb.obs,
            design_factors=[groupby, "gender", "manual_anno_L2"])
        sc.pp.filter_genes(dds, min_cells=1)
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=[groupby, vals[0], vals[1]])
        stat_res.summary()
        res = stat_res.results_df
        res['Symbol'] = res.index
        vres = res[(abs(res.log2FoldChange) < 2.5)] # Remove some weird outlier
        volcano(vres, symbol='Symbol',
                save=cmf.ensure_file(f"{fig_folder}{cluster}/{groupby}"),
                colors=['red', 'lightgrey', 'red'],
                log2fc_thresh=0.5
                )
        res.to_csv(cmf.ensure_file(f"{fig_folder}{cluster}/deseq_results.tsv"), sep="\t")
        # GSEA
        ranking = res[['Symbol', 'stat']].dropna().sort_values('stat', ascending = False)
        ranking = ranking.drop_duplicates('Symbol')
        ranking.index = ranking.index.str.upper()
        for gs in params.gene_sets:
            pre_res = gseapy.prerank(rnk=ranking.drop(columns=['Symbol']),
                                 gene_sets=gs, threads=24)
            out = []
            for term in list(pre_res.results):
                out.append([term,
                            pre_res.results[term]['fdr'],
                            pre_res.results[term]['es'],
                            pre_res.results[term]['nes']])
            out_df = pd.DataFrame(out, columns=['term', 'fdr', 'es', 'nes']).sort_values('nes', ascending=False).reset_index(drop=True)
            out_df = out_df[out_df['fdr'] < 0.5]
            out_df = out_df.head(30)
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
            plt.savefig(cmf.ensure_file(f"{fig_folder}{cluster}/GSEA_genesets/{gs}.png"), dpi=300)
            plt.close(fig)
            plt.clf()
            for term in pre_res.results.keys():
                if "ALS" in term or "AMYOTROPHIC" in term:
                    gseaplot(figsize=(6, 10), ofname=cmf.ensure_file(f"{fig_folder}{cluster}/GSEA_ALS/plot_{groupby}_{gs}_{term}.png"),
                             term=term,
                             **pre_res.results[term])

        pre_res = gseapy.prerank(rnk=ranking.drop(columns=['Symbol']),
                                 gene_sets="KEGG_2021_Human")
        pre_res.res2d['Obs'] = cluster
        enr_results.append(pre_res.res2d)

        sigs = res[res.padj < 0.05].copy()
        dedf = sigs.copy().sort_values('log2FoldChange', ascending=False).reset_index(drop=True)

        sig_size[cluster] = len(sigs)
        sigs['log2FoldChange'] = np.abs(sigs['log2FoldChange'])
        sigs.sort_values('log2FoldChange', ascending=False, inplace=True)
        genes_heatmap[cluster] = sigs.index.tolist()[:10]

        # Only female from now on!
        adata_c = adata_c[adata_c.obs["gender"] == "Female"].copy()
        adata_small = adata_c[:, sigs.index.tolist()[:100]]
        sc.pp.normalize_total(adata_small)
        sc.pp.log1p(adata_small)
        sc.pp.scale(adata_small)
        sc.pp.neighbors(adata_small, n_neighbors=10, use_rep="X")
        sc.tl.umap(adata_small)
        sc.pl.umap(adata_small, color=groupby, show=False)
        plt.savefig(cmf.ensure_file(f"{fig_folder}{cluster}/DE_Female_only.png"), dpi=300)
        top_gene = sigs.index.tolist()[0]
        sc.pl.umap(adata_small, color=top_gene, show=False)
        plt.savefig(cmf.ensure_file(f"{fig_folder}{cluster}/DE_{top_gene}_Female_only.png"), dpi=300)
        # sc.pl.umap(adata_small, color="MYC", show=False)
        # plt.savefig(cmf.ensure_file(f"{fig_folder}{cluster}/DE_MYC_Female_only.png"), dpi=300)

        # Top and bottom ton n log fold genes heatmaps
        top_n = 10

        adata_h = adata_c.copy()
        sc.pp.normalize_total(adata_h)
        sc.pp.log1p(adata_h)
        sc.pp.scale(adata_h)
        for gene in genes_heatmap[cluster]:
            i = np.where(adata_h.var_names == gene)[0][0]
            a = adata_h[adata_h.obs[groupby] == vals[0]].X[:,i].todense()
            b = adata_h[adata_h.obs[groupby] == vals[1]].X[:,i].todense()
            p_value = stats.mannwhitneyu(np.asarray(a), np.asarray(b)).pvalue
            fig, ax = plt.subplots(figsize=(7, 5))
            sc.pl.violin(adata_h, gene, groupby=groupby, ax=ax, show=False, use_raw=False)
            ax.text(x=0.5, y=max(ax.get_ylim()),
                    s=f'p-value: {p_value}', ha='center', va='bottom')
            plt.savefig(cmf.ensure_file(f"{fig_folder}{cluster}/genes/{cluster}_{gene}.png"), dpi=300)
            plt.close(fig)
            plt.clf()
        genes_to_show = dedf[-top_n:].Symbol.tolist() + dedf[:top_n].Symbol.tolist()
        sc.pl.heatmap(adata_h, genes_to_show, groupby=groupby, swap_axes=True, show=False, use_raw=False)
        plt.savefig(cmf.ensure_file(f"{fig_folder}{cluster}/top{top_n}_de_genes.png"), dpi=300)
        if iter_obs == params.cell_type:
            adata_h = adata_c.copy()
            pbs = []
            for i in adata_h.obs["individual"].unique():
                samp_cell_subset = adata_h[adata_h.obs["individual"] == i]
                rep_adata = sc.AnnData(X=samp_cell_subset.X.sum(axis=0),
                                       var=samp_cell_subset.var[[]])
                rep_adata.obs_names = [i]
                rep_adata.obs[groupby] = samp_cell_subset.obs[groupby].iloc[0]
                pbs.append(rep_adata)
            pb = sc.concat(pbs)
            pb.X = np.array(pb.X)
            sc.pp.normalize_total(pb)
            sc.pp.log1p(pb)
            sc.pp.scale(pb)
            sc.pl.heatmap(pb, genes_to_show, groupby=groupby, swap_axes=True, show=False, use_raw=False)
            plt.savefig(cmf.ensure_file(f"{fig_folder}{cluster}/top{top_n}_de_genes_pb.png"), dpi=300)
    except Exception as e:
        traceback.print_exc()
        print(f"{cluster} GSEA not performed")

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata)

# Concatenate the GSEA results for each cluster into a single DataFrame
enr_combined = pd.concat(enr_results)

# Pivot the DataFrame to prepare for heatmap plotting
heatmap_data = enr_combined.pivot(index='Term', columns='Obs', values='NES')
heatmap_data = heatmap_data.fillna(0)
top_pathways = heatmap_data.max(axis=1).nlargest(100).index
heatmap_data = heatmap_data.loc[top_pathways]
heatmap = sns.clustermap(
    heatmap_data,
    cmap='RdYlBu_r',
    cbar_kws={'label': 'NES'},
    yticklabels=True,
    figsize=(16, 24)
)

# Set labels and title
heatmap.ax_heatmap.set_xlabel('Cluster')
heatmap.ax_heatmap.set_ylabel('Pathway')
heatmap.ax_heatmap.set_title('GSEA Heatmap ')

# Rotate and align the X-axis labels
plt.setp(heatmap.ax_heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
plt.tight_layout()
# Save the figure to a PNG file
heatmap.savefig(cmf.ensure_file(f'{fig_folder}Celltype_{groupby}_GSEA_heatmap.png'), dpi=300)

# Significant DEG number barplot
plt.clf()
df = pd.DataFrame(list(sig_size.items()), columns=['Obs', 'Number'])
sns.barplot(x='Obs', y='Number', data=df)
plt.title('Number of significant DEG per Obs')
plt.savefig(f'{fig_folder}sig_barplot.png', dpi=300)
plt.clf()

# DEG heatmap
plt.clf()
sc.pl.heatmap(adata, genes_heatmap, groupby=groupby, dendrogram=True, use_raw=False,
                  show_gene_labels=True, swap_axes=True, show=False)
plt.savefig(cmf.ensure_file(f'{fig_folder}de_heatmap.png'), dpi=300)
plt.clf()


