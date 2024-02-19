import anndata as ad
import pandas as pd
from params import Params
import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.plots import volcano
import numpy as np
import matplotlib

matplotlib.use('agg')

params = Params()

adata = ad.read_h5ad(params.folder + "tcre_adata.h5ad")
print("Loaded adata")

groupby = ["treatment", "status"][1]
np.unique(adata.obs[groupby])
if groupby == "treatment":
    vals = ["Ctrl", "Ropi"]
    print(len(adata))
    adata = adata[adata.obs["status"] == "ALS"].copy()
    adata = adata[adata.obs["day"].isin(["D7", "D14"])].copy()
    print(len(adata))
elif groupby == "status":
    vals = ["healthy", "ALS"]
    print(len(adata))
    # adata = adata[adata.obs["treatment"] != "Ropi"].copy()
    print(len(adata))

pbs = []
column_name_level = "manual_anno_L2"
for individual in adata.obs["individual"].unique():
    adata_i = adata[adata.obs['individual'] == individual]
    for gr in adata_i.obs[groupby].unique():
        adata_grb = adata_i[adata_i.obs[groupby] == gr]
        for cluster in adata_grb.obs[column_name_level].unique():
            samp_cell_subset = adata_grb[adata_grb.obs[column_name_level] == cluster]

            rep_adata = sc.AnnData(X=samp_cell_subset.X.sum(axis=0),
                                   var=samp_cell_subset.var[[]])
            rep_adata.obs_names = [individual + "_" + cluster + "_" + gr]
            rep_adata.obs[groupby] = samp_cell_subset.obs[groupby].iloc[0]
            rep_adata.obs['gender'] = samp_cell_subset.obs['gender'].iloc[0]
            rep_adata.obs[column_name_level] = samp_cell_subset.obs[column_name_level].iloc[0]
            pbs.append(rep_adata)
pb = sc.concat(pbs)

# adata_dense = np.asarray(pb.X)
# condition_1_mask = np.sum(adata_dense >= 5, axis=0) >= 10
# total_counts_per_sample = np.sum(adata_dense, axis=1)
# cpm = (adata_dense / total_counts_per_sample[:, np.newaxis]) * 1e6
# condition_2_mask = np.sum(cpm >= 0.5, axis=0) >= 10
# genes_to_keep = condition_1_mask & condition_2_mask
# counts = pd.DataFrame(adata_dense[:, genes_to_keep], columns=np.array(pb.var_names)[genes_to_keep])
counts = pd.DataFrame(pb.X, columns=pb.var_names)

dds = DeseqDataSet(
    counts=counts,
    metadata=pb.obs,
    design_factors=[groupby, "gender", column_name_level])
sc.pp.filter_genes(dds, min_cells=1)
dds.deseq2()
# sc.tl.pca(dds)
# sc.pl.pca(dds, color='status', size=200, show=True, save="deseq2.png")
stat_res = DeseqStats(dds, contrast=[groupby, vals[0], vals[1]])
stat_res.summary()
res = stat_res.results_df
res['Symbol'] = res.index
res.to_csv("deseq_results.tsv", sep="\t")

vres = res[(abs(res.log2FoldChange) < 2.5)]
volcano(vres, symbol='Symbol',
        save=f"deseq_volcano",
        colors=['red', 'lightgrey', 'red'],
        log2fc_thresh=0.5
        )
