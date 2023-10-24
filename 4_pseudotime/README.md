# Braak-stage-informed pseudotime

To reproduce the core results (pseudotime) used for results shown in the paper, the following
order of scripts can be run

1. `preprocessing.py`: Reduces file size by removing information redundant for pseudotime
    construction and CellRank-based analyses
2. `tmaps.py`: Computes transport maps with moscot
3. `dpt.py`: Computes diffusion pseudotime
4. `cr.py`: Runs CellRank pipeline to identify terminal states and driver genes

Note: we provide examples for `adapt_homeo_analysis` here. The same analysis can be applied to the ADAM and PVM subsets by replacing `MG_Adapt` and `MG_Homeo` with `MG_ADAM` and `MG_PVM`, respectively.