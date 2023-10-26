import scanpy as sc
import pandas as pd
import numpy as np

# Load datasets
fmg = sc.read_h5ad('/path/to/FreshMG.h5ad')
pad = sc.read_h5ad('/path/to/PsychAD.h5ad')

# Filter to AD and CTRL samples
fmg = fmg[fmg.obs['dx'].isin(['AD','CTRL'])]  
pad = pad[pad.obs['dx'].isin(['AD','CTRL'])]

# Clean datasets
def clean_adata(adata):
    adata_clean = sc.AnnData(adata.X, adata.obs_names.to_frame(), adata.var_names.to_frame())
    adata_clean.obs['batch'] = adata.obs['batch'] 
    return adata_clean

fmg = clean_adata(fmg)
pad = clean_adata(pad)

# Concatenate datasets
combined = ad.concat({"FreshMG": fmg, "PsychAD": pad}, label="batch") 

# Identify highly variable genes
sc.pp.normalize_total(combined, target_sum=1e4)
sc.pp.log1p(combined)
sc.pp.highly_variable_genes(combined, min_mean=0.0125, max_mean=3, batch_key='batch')

# Filter to HVG
combined = combined[:, combined.var['highly_variable']]
fmg = fmg[:, combined.var_names]
pad = pad[:, combined.var_names]

# Write output files
fmg.write_h5ad('freshmg_clean.h5ad')
pad.write_h5ad('psychad_clean.h5ad')
combined.write_h5ad('combined_clean.h5ad')
combined.to_df().to_csv('combined_clean.csv')