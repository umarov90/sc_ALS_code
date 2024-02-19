import anndata as ad
import pandas as pd
import gc
import psutil
from params import Params
import dask.dataframe as dd


def get_human_readable(size, precision=2):
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB']
    suffixIndex = 0
    while size > 1024 and suffixIndex < 4:
        suffixIndex += 1  # increment the index of the suffix
        size = size / 1024.0  # apply the division
    return "%.*f%s" % (precision, size, suffixes[suffixIndex])


def print_memory():
    mem = psutil.virtual_memory()
    print(f"used: {get_human_readable(mem.used)} available: {get_human_readable(mem.available)}")


if __name__ == '__main__':
    params = Params()
    adata = ad.read_h5ad(params.folder + "tcre_adata.h5ad")
    all_dfs = {}
    column_name_level = "manual_anno_L2"
    for sample_id in adata.obs["sample_id"].unique():
        print(sample_id)
        adata_s = adata[adata.obs["sample_id"] == sample_id]
        ctss_path = f"/osc-fs_home/hon-chun/analysis/ramzan/ALS/scafe/solo/{sample_id}/bam_to_ctss/{sample_id}/bed/{sample_id}.CB.ctss.bed.gz"
        sm = pd.read_csv(ctss_path, sep="\t", header=None, usecols=[0, 1, 2, 3, 4])
        sm.columns = ["chr", "start", "end", "barcode", "count"]
        print(f"Max barcode: {sm['barcode'].str.len().max()}")
        print(len(sm))

        gc.collect()
        print_memory()
        
        # Step 1: Create the barcode-cluster mapping DataFrame
        barcode_cluster_list = []
        clusters = adata_s.obs[column_name_level].unique()
        for cluster in clusters:
            adata_c = adata_s[adata_s.obs[column_name_level] == cluster]
            barcodes = ('\t'.join(adata_c.obs['old_obs_names'])).split('\t')
            barcodes = [x + "-1" for x in barcodes]
            barcodes = [x.split('_', 1)[-1] for x in barcodes]
            for barcode in barcodes:
                barcode_cluster_list.append({'barcode': barcode, 'cluster': cluster})

        barcode_cluster_df = pd.DataFrame(barcode_cluster_list)

        # Convert 'barcode_cluster_df' to Dask DataFrame for efficient merging
        barcode_cluster_dask_df = dd.from_pandas(barcode_cluster_df, npartitions=24)

        # Step 2: Merge with 'sm' DataFrame
        dask_df = dd.from_pandas(sm, npartitions=24)
        merged_df = dask_df.merge(barcode_cluster_dask_df, on='barcode', how='inner')

        # Compute the merged DataFrame (optional step, depending on subsequent needs)
        merged_df = merged_df.compute()

        # Step 3: Split the merged DataFrame by cluster and update 'all_dfs'
        for cluster in clusters:
            cluster_df = merged_df[merged_df['cluster'] == cluster].copy()
            cluster_df.drop(['barcode', 'cluster'], axis=1, inplace=True)  # Optional: Drop the 'barcode' column if not needed
            all_dfs.setdefault(cluster, []).append(cluster_df)


        # for cluster in adata_s.obs[column_name_level].unique():
        #     print(cluster)
        #     adata_c = adata_s[adata_s.obs[column_name_level] == cluster]
        #     barcodes = ('\t'.join(adata_c.obs['old_obs_names'])).split('\t')
        #     barcodes = [x + "-1" for x in barcodes]
        #     barcodes = [x.split('_', 1)[-1] for x in barcodes]
        #     print(len(barcodes))
        #     dask_df = dd.from_pandas(sm, npartitions=32)
        #     filtered_dask_df = dask_df[dask_df['barcode'].isin(barcodes)]
        #     filtered_df = filtered_dask_df.compute()
        #     filtered_df.drop('barcode', axis=1, inplace=True)
        #     print(len(filtered_df))
        #     all_dfs.setdefault(cluster, []).append(filtered_df)

    for key in all_dfs.keys():
        concatenated_df = pd.concat(all_dfs[key], axis=0, ignore_index=True)
        concatenated_df.to_csv(params.folder + "tracks/" + key + ".bed.gz",
                               compression="gzip", sep="\t", index=False, header=False)
        print(f"Saved {key}")
        