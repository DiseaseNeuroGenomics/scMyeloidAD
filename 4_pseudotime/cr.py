from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import load_npz

import cellrank as cr
import scanpy as sc


# General settings
cr.settings.verbosity = 2
sc.settings.verbosity = 2

# Constants
DATA_PATH = Path("../data/")

STAGE_KEY = 'braak'

CELLTYPE_SUBSET = ["MG_Adapt", "MG_Homeo"]


if __name__ == "__main__":
    adata = sc.read(DATA_PATH / "processed" / "adata.h5ad")
    celltype_palette = dict(zip(adata.obs["subclass"].cat.categories, adata.uns["subclass_colors"]))
    celltype_fine_palette = dict(zip(adata.obs["subtype"].cat.categories, adata.uns["subtype_colors"]))

    adata = adata[adata.obs["subclass"].isin(CELLTYPE_SUBSET), :].copy()

    dpt_df = pd.read_csv(
        DATA_PATH / "results" / "dpt" / "adapt_homeo.csv",
        index_col=0
    )

    percentile = 95
    threshold = np.percentile(dpt_df["dpt_pseudotime"], q=percentile)
    obs_mask = adata.obs_names.isin(dpt_df.loc[dpt_df["dpt_pseudotime"] < threshold].index)
    dpt_df = dpt_df.loc[adata.obs_names[obs_mask]]

    adata = adata[obs_mask, :].copy()
    adata.obs["dpt_pseudotime"] = dpt_df["dpt_pseudotime"].values

    adata.obs[f"subclass_{STAGE_KEY}"] = adata.obs["subclass"].astype(str) + "-" + adata.obs[STAGE_KEY].astype(str)
    adata.obs[f"subclass_{STAGE_KEY}"] = adata.obs[f"subclass_{STAGE_KEY}"].astype("category")

    adata.obs[f"subtype_{STAGE_KEY}"] = adata.obs["subtype"].astype(str) + "-" + adata.obs[STAGE_KEY].astype(str)
    adata.obs[f"subtype_{STAGE_KEY}"] = adata.obs[f"subtype_{STAGE_KEY}"].astype("category")

    adata.uns["subclass_colors"] = celltype_palette
    adata.uns["subtype_colors"] = celltype_fine_palette

    # CellRank analysis
    stages = np.sort(adata.obs[STAGE_KEY].unique())

    wk = cr.external.kernels.WOTKernel(adata, time_key=STAGE_KEY)

    tmat = load_npz(DATA_PATH / "results" / "cr" / "tmat_adapt_homeo.npz")
    tmat = tmat[obs_mask, :]
    tmat = tmat[:, obs_mask].copy()
    wk.transition_matrix = tmat

    estimator = cr.estimators.GPCCA(wk)
    estimator.compute_schur(n_components=5)

    adata.uns["subclass_colors"] = [adata.uns["subclass_colors"][celltype] for celltype in adata.obs["subclass"].cat.categories]
    adata.uns["subtype_colors"] = [adata.uns["subtype_colors"][celltype] for celltype in adata.obs["subtype"].cat.categories]

    estimator.compute_macrostates(n_states=5, cluster_key="subclass")
    estimator.set_terminal_states(states=["MG_Adapt_1, MG_Adapt_2, MG_Adapt_3, MG_Adapt_4", "MG_Homeo"])
    estimator.rename_terminal_states({"MG_Adapt_1, MG_Adapt_2, MG_Adapt_3, MG_Adapt_4": "MG_Adapt"})

    estimator.compute_absorption_probabilities()
    drivers = estimator.compute_lineage_drivers(
        lineages=["MG_Adapt", "MG_Homeo"],
        cluster_key=STAGE_KEY,
        clusters=[0.0],
        return_drivers=True,
    )

    drivers.to_csv(DATA_PATH / "results" / "cr" / "drivers_adapt_homeo.csv")