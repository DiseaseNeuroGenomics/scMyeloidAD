import os
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz, save_npz

import cellrank as cr
import scanpy as sc
import scvelo as scv
from scanpy.tools._dpt import DPT


# General setting
cr.settings.verbosity = 2
sc.settings.verbosity = 2

# Constants
DATA_PATH = Path("../data/")

STAGE_KEY = 'braak'

CELLTYPE_SUBSET = ["MG_Adapt", "MG_Homeo"]

# Function definitions
def get_symmetric_transition_matrix(transition_matrix):
    sym_mat = (transition_matrix + transition_matrix.T) / 2

    # normalise transition matrix
    row_sums = sym_mat.sum(axis=1).A1
    sym_mat.data = sym_mat.data / row_sums[sym_mat.nonzero()[0]]

    return sym_mat


def set_transition_matrix(wot_kernel, tmap, last_time_point, conn_weight, threshold="auto"):
    if threshold:
        wot_kernel._threshold_transport_maps(tmap, threshold)
    tmap = wot_kernel._restich_tmaps(tmap, last_time_point, conn_weight=conn_weight)
    wot_kernel.transition_matrix = tmap.X


if __name__ == "__main__":
    os.makedirs(DATA_PATH / "results" / "dpt", exist_ok=True)

    adata = sc.read(DATA_PATH / "processed" / "adata.h5ad")

    adata = adata[adata.obs["subclass"].isin(CELLTYPE_SUBSET), :].copy()
    
    # CellRank analysis
    stages = np.sort(adata.obs[f"{STAGE_KEY}"].unique())

    wk = cr.external.kernels.WOTKernel(adata, time_key=f"{STAGE_KEY}")
    if (DATA_PATH / "results" / "cr" / "tmat_adapt_homeo.npz").is_file():
        wk.transition_matrix = load_npz(DATA_PATH / "results" / "cr" / "tmat_adapt_homeo.npz")
    else:
        transport_maps = {
            key: sc.read(DATA_PATH / "results" / "tmaps" / f"tmap_{key}_{'-'.join(CELLTYPE_SUBSET)}.h5ad")
            for key in zip(stages[:-1], stages[1:])
        }

        set_transition_matrix(
            wot_kernel=wk, tmap=transport_maps, last_time_point="all", conn_weight=0.1
        )

        del transport_maps
        save_npz(
            DATA_PATH / "results" / "cr" / "tmat_adapt_homeo.npz", wk.transition_matrix
        )

    # Pseudotime construction
    dpt = DPT(adata=adata, neighbors_key='neighbors')
    dpt._transitions_sym = get_symmetric_transition_matrix(wk.transition_matrix)
    dpt.compute_eigen(n_comps=15, random_state=0)

    adata.obsm['X_diffmap'] = dpt.eigen_basis
    adata.uns['diffmap_evals'] = dpt.eigen_values

    adata.uns['iroot'] = 236768  # See dpt_root_id.ipynb
    sc.tl.dpt(adata)

    adata = adata[~(adata.obs["dpt_pseudotime"] == np.inf), :].copy()

    adata.obs[["dpt_pseudotime", STAGE_KEY, "subclass", "subtype"]].to_csv(
        DATA_PATH / "results" / "dpt" / "adapt_homeo.csv"
    )

    percentile = 95
    threshold = np.percentile(dpt["dpt_pseudotime"], q=percentile)
    adata = adata[adata.obs_names[adata.obs["dpt_pseudotime"] < threshold], :].copy()
    adata.obs[["dpt_pseudotime", STAGE_KEY, "subclass", "subtype"]].to_csv(
        DATA_PATH / "results" / "dpt" / "adapt_homeo-outliers_removed.csv"
    )
