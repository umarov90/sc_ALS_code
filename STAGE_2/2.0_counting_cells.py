import os
import sys
import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
from scipy import stats
from gseapy.plot import gseaplot
from sanbomics.plots import volcano
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
matplotlib.use('agg')
script_dir = os.path.dirname(__file__)  # Directory containing myscript.py
parent_dir = os.path.dirname(script_dir)  # Project directory
sys.path.append(parent_dir)
from utils.params import Params

params = Params()
adata = ad.read_h5ad(params.file_path)
adata = adata.raw.to_adata()

num_tot_cells = adata.obs.groupby(['individual']).count()
num_tot_cells = dict(zip(num_tot_cells.index, num_tot_cells.day))
# does not matter if day or pool, all are the same
cell_type_counts = adata.obs.groupby(['individual', 'status', 'manual_anno_L0']).count()
cell_type_counts = cell_type_counts[cell_type_counts.sum(axis = 1) > 0].reset_index()
cell_type_counts['total_cells'] = cell_type_counts.individual.map(num_tot_cells).astype(int)
cell_type_counts['frequency'] = cell_type_counts.day / cell_type_counts.total_cells

p_values = {}
cell_types = cell_type_counts['manual_anno_L0'].unique()
for cell_type in cell_types:
    subset = cell_type_counts[cell_type_counts['manual_anno_L0'] == cell_type]
    freq_healthy = subset[subset['status'] == 'healthy']['frequency']
    freq_ALS = subset[subset['status'] == 'ALS']['frequency']
    p_value = stats.mannwhitneyu(np.asarray(freq_healthy), np.asarray(freq_ALS)).pvalue
    p_values[cell_type] = p_value

plt.figure(figsize=(10, 4))
ax = sns.boxplot(data=cell_type_counts, x='manual_anno_L0', y='frequency', hue='status',
                 hue_order=["healthy", "ALS"])
plt.xticks(rotation=35, rotation_mode='anchor', ha='right')

# Adjust this factor to move the annotation up or down
offset_factor = 0.05
# Get the current y-axis limits to help position the p-value annotations dynamically
y_min, y_max = ax.get_ylim()
offset = (y_max - y_min) * offset_factor

# Annotate plot with p-values
for i, cell_type in enumerate(cell_types):
    y_pos = y_max + offset
    plt.text(i, y_pos, f'p={p_values[cell_type]:.2e}', horizontalalignment='center', size='small')

plt.savefig(f'counts.png', dpi=300)
