# Python script to extract reduce AnnData object to information relevant for pseudotime-based analysis
import os
from pathlib import Path

import scanpy as sc

# Constants
# Path to directory containing data
DATA_PATH = Path("./data/")

# Name of H5AD file
fname = "data.h5ad"

# Subset of interest for pseudotime-based analysis
CELLTYPE_SUBSET = ["MG_Adapt", "MG_Homeo", "MG_PVM", "MG_ADAM"]
STAGE_KEY = "braak"

OBS_COLS_TO_KEEP = ["Donor", "age", "sex", "race", "class", "subclass", "subtype", "CERAD", "Braak", "dementia"]
VAR_COLS_TO_KEEP = ["highly_variable_features"]


if __name__ == "__main__":
    os.makedirs(DATA_PATH / "processed", exist_ok=True)

    # Data loading
    adata = sc.read(DATA_PATH / fname)

    cluster_to_color = dict(zip(adata.obs["subclass"].cat.categories, adata.uns["subclass_colors"]))
    subtype_to_color = dict(zip(adata.obs["subtype"].cat.categories, adata.uns["subtype_colors"]))

    # Data cleaning
    adata.obs.drop(adata.obs.columns.difference(OBS_COLS_TO_KEEP), axis=1, inplace=True)
    adata.obs.columns = adata.obs.columns.str.lower()

    adata.var.drop(adata.var.columns.difference(VAR_COLS_TO_KEEP), axis=1, inplace=True)
    adata.var.rename(columns={"highly_variable_features": "highly_variable"}, inplace=True)

    adata.uns = {}
    adata.uns["subclass_colors"] = [cluster_to_color[cluster] for cluster in adata.obs["subclass"].cat.categories]
    adata.uns["subtype_colors"] = [subtype_to_color[cluster] for cluster in adata.obs["subtype"].cat.categories]

    for key in [
        "X_pca",
        "X_pca_regressed",
        "X_umap",
        "pca_regressed_harmony_knn_distances",
        "pca_regressed_harmony_knn_indices",
    ]:
        del adata.obsm[key]

    adata.obsm["X_pca"] = adata.obsm.pop("X_pca_regressed_harmony")

    del adata.varm["de_res"]
    del adata.obsp["W_pca_regressed_harmony"]

    adata = adata[adata.obs["class"].isin(CELLTYPE_SUBSET), adata.var["highly_variable"]]
    adata = adata[adata.obs[STAGE_KEY].notnull(), :].copy()

    # Data preprocessing
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=15, use_rep="X_pca")

    # Data writing
    adata.write(DATA_PATH / "processed" / "adata.h5ad")
