import scanpy as sc
import pandas as pd
import numpy as np

RESULTS_FOLDER = '/path/to/results' 

ANNDATA_FNAME = os.path.join(RESULTS_FOLDER, 'anndata.h5ad')
CSV_FNAME = os.path.join(RESULTS_FOLDER, 'combined.csv')

ADJ_FNAME = os.path.join(RESULTS_FOLDER, 'adjacencies.tsv')  
MODULES_FNAME = os.path.join(RESULTS_FOLDER, 'modules.csv')
REGULONS_FNAME = os.path.join(RESULTS_FOLDER, 'regulons.csv')

#Reference pyscenic documentation for latest database
db_fpath = 'path/to/databases'
f_db = os.path.join(db_fpath, 'genes_vs_motifs.feather')
motif_annotations = os.path.join(db_fpath, 'motif_annotations.tbl')
tf_list = os.path.join(db_fpath, 'tf_list.txt')


# File paths from preprocessing

FRESHMG_FNAME = 'freshmg_clean.h5ad' 
PSYCHAD_FNAME = 'psychad_clean.h5ad'

GRNBOOST_OUT = 'grnboost_out.csv'
PYSCENIC_ADJ = 'pyscenic_adj.tsv'
PYSCENIC_REGULONS = 'pyscenic_regulons.csv'
PYSCENIC_AUCELL = 'pyscenic_aucell.csv'

# Prepare input data 

df = pd.concat([
   sc.read_h5ad(FRESHMG_FNAME).to_df(),
   sc.read_h5ad(PSYCHAD_FNAME).to_df()
])

df.to_csv('grn_input.csv')

# Run GRNBoost2
!arboreto_with_multiprocessing.py $CSV_FNAME $tf_list 
   --method grnboost2  
   --output $ADJ_FNAME
   --num_workers 20
   --seed 123

# Run pyscenic ctx
!pyscenic ctx $ADJ_FNAME $f_db  
   --annotations_fname $motif_annotations
   --expression_mtx_fname $CSV_FNAME 
   --output $MODULES_FNAME
   
# Run pyscenic aucell   
!pyscenic aucell $CSV_FNAME $MODULES_FNAME
   --output $REGULONS_FNAME