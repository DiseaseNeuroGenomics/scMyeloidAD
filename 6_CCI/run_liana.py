#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Infering cell-cell interactions with LIANA
"""
import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
import logging
import traceback

adata = sc.read_h5ad(data_path, backed="r")

for i, dnr in enumerate(donors):
    adata_sub = adata[adata.obs[sample_label] == donors[i]].to_memory()

    try:
        li.mt.rank_aggregate(
            adata_sub,
            groupby=cell_label,
            resource_name="consensus",
            expr_prop=0.1,  # must be expressed in expr_prop fraction of cells
            use_raw=False,  # run on log- and library-normalized counts
            min_cells=5,
            n_perms=1000,
            verbose=True,
            inplace=True,
        )

        adata_sub.uns["liana_res"].to_csv(results_path, index=False)
    except Exception as e:
        logging.error(traceback.format_exc())
        continue
