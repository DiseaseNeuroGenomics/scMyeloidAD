# sc
import pegasus as pg
import scanpy as sc
import anndata as ad
import pge

# plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

# data
import numpy as np
import pandas as pd
from scipy import stats
from scipy import sparse
import h5py

### parameters
csv_manifest = '/path/to/manifest.csv'
list_attributes = ['Donor','Source','Type','Inhibitors','Sequencing']
genome = 'GRCh38'
n_genes_lower = 500
n_counts_lower = 1000
n_genes_upper = 8000
n_counts_upper = 40000
percent_cells=0.05
percent_mito=20

### aggregate
data = pg.aggregate_matrices(csv_file=csv_manifest, attributes=list_attributes, default_ref=genome, append_sample_name=True)

### identify robust genes
pg.identify_robust_genes(data, percent_cells=percent_cells)

### remove features that are not robust (expressed at least 0.05% of cells) and protein_coding from downstream analysis
data._inplace_subset_var(data.var['robust'])

### nUMI + nGene based QC
pg.qc_metrics(data, min_genes=n_genes_lower, max_genes=n_genes_upper, min_umis=n_counts_lower, max_umis=n_counts_upper, mito_prefix='MT-', percent_mito=percent_mito)
pg.filter_data(data)

### clean unused categories
data = pge.clean_unused_categories(data)

### log normalization
pg.log_norm(data)