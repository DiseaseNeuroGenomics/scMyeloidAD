import os
from pathlib import Path, PosixPath

from tqdm import tqdm

import numpy as np
import pandas as pd

import scanpy as sc
from anndata import AnnData
from moscot.problems import TemporalProblem


# General settings
sc.settings.verbosity = 2


# Function definitions
def save_transport_maps(time_points, problem, adata: AnnData, stage_key: str, fpath: PosixPath, fname_appendix: str):
    """Save transport maps computed with moscot."""
    for t0, t1 in tqdm(zip(time_points[:-1], time_points[1:])):
        AnnData(
            X=np.array(problem.solutions[(t0, t1)].to("cpu").transport_matrix),
            obs=pd.DataFrame(index=adata.obs_names[adata.obs[stage_key] == t0]),
            var=pd.DataFrame(index=adata.obs_names[adata.obs[stage_key] == t1]),
        ).write(fpath / f"tmap_{(t0, t1)}_{'-'.join(fname_appendix)}.h5ad")


# Constants
# Path to directory containing data
DATA_PATH = Path("../data/")

STAGE_KEY = "braak"
ORGANISM = "human"

CELLTYPE_SUBSET = ["MG_Adapt", "MG_Homeo"]

# Whether or not to save transport maps
SAVE = True


if __name__ == "__main__":
    os.makedirs(DATA_PATH / "results" / "tmaps", exist_ok=True)

    # Data loading
    adata = sc.read(DATA_PATH / "processed" / "adata.h5ad")
    adata = adata[adata.obs["subclass"].isin(CELLTYPE_SUBSET), :].copy()

    # Data preprocessing
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=15, use_rep="X_pca")

    # Optimal transport analysis
    epsilon = 0.05

    tau_a = 1
    tau_b = 1

    tp = TemporalProblem(adata)
    tp.score_genes_for_marginals()
    tp.prepare(f"{STAGE_KEY}", joint_attr="X_pca")

    stages = np.sort(adata.obs[f"{STAGE_KEY}"].unique())


    tp.solve(
        epsilon=epsilon,
        tau_a=tau_a,
        tau_b=tau_b,
        batch_size=1024,
        scale_cost="mean",
        device="gpu",
    )

    if SAVE:
        save_transport_maps(
            time_points=stages,
            problem=tp,
            adata=adata,
            stage_key=STAGE_KEY,
            fpath=DATA_PATH / "tmaps",
            fname_appendix=CELLTYPE_SUBSET,
        )
