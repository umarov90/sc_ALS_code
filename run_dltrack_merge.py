import anndata as ad
import sys
import pandas as pd
import gc
import psutil
import glob
import os
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
    # Step 1: List all .bed.gz files
    files = glob.glob(params.folder + "tracks/" + '*.bed.gz')
    # Step 2: Create a dictionary to hold data for each group
    data_groups = {}
    for file in files:
        print(file)
        # Extract the group identifier (between the first "_" and ".bed.gz")
        basename = os.path.basename(file)
        group_id = basename.split('_', 1)[1].rsplit('.bed.gz', 1)[0]
        print(group_id)
        # Read the .bed.gz file
        df = pd.read_csv(file, sep='\t', compression='gzip', header=None)
        data_groups.setdefault(group_id, []).append(df)

    # Step 3: Concatenate and save each group's DataFrames
    for group_id, dfs in data_groups.items():
        print(group_id)
        # Concatenate all DataFrames in the list
        concatenated_df = pd.concat(dfs, ignore_index=True)
        
        # Define the output filename and path
        output_path = f"{params.folder}dl_tracks/{group_id}.bed.gz"
        
        # Save the concatenated DataFrame to .bed.gz format
        concatenated_df.to_csv(output_path, sep='\t', index=False, header=False, compression='gzip')

    print("Files have been processed and saved.")