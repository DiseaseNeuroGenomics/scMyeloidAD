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
prefix='example_run'
flavor='cell_ranger'
batch_key='Source'
n_top_genes=1000
npc=50
K=100
n_neighbors=50
res=2.0

### load QC norm data
data = pg.read_input(prefix+'_qc_norm.h5ad')

### remove features that are not robust from downstream analysis
data._inplace_subset_var(data.var['robust'])

### HVF
pge.scanpy_hvf(data, flavor=flavor, batch_key=batch_key, n_top_genes=n_top_genes, robust_protein_coding=True)

### pca/harmony/umap/tsne
pg.pca(data, n_components=npc)
pg.elbowplot(data)
npc = data.uns["pca_ncomps"]
pg.regress_out(data, attrs=['n_counts','percent_mito','cycle_diff'])
pg.run_harmony(data, batch=batch_key, rep='pca_regressed', max_iter_harmony=20, n_comps=npc)
pg.neighbors(data, rep='pca_regressed_harmony', use_cache=False, dist='cosine', K=K, n_comps=npc)
pg.umap(data, rep='pca_regressed_harmony', n_neighbors=n_neighbors, rep_ncomps=npc)
pg.tsne(data, rep='pca_regressed_harmony', rep_ncomps=npc)

### clustering
pg.leiden(data, rep='pca_regressed_harmony', resolution=res, class_label='leiden_labels_res20')

# de
pg.de_analysis(data, cluster='leiden_labels_res20')

# marker
marker_dict_res20 = pg.markers(data)
pg.write_results_to_excel(marker_dict_res20, prefix+"_pass2_res20.de.xlsx")

# annotate
celltype_dict_res20 = pg.infer_cell_types(data, markers = 'human_brain')
cluster_names_res20 = pg.infer_cluster_names(celltype_dict_res20)
pg.annotate(data, name='anno_res20', based_on='leiden_labels_res20', anno_dict=cluster_names_res20)

### plot umap
pg.scatter(data, attrs=['leiden_labels_res20','anno_res20'], basis='umap', legend_loc='on data', wspace=0.1)