# scMyeloidAD
Plasticity of Human Microglia and Brain Perivascular Macrophages in Aging and Alzheimer’s Disease

# Contents

```
scMyeloidAD
├──1_taxonomy: preprocessing and QC of single-cell data
├──2_crumblr: analysis on compositional variation
├──3_dreamlet: analysis on gene expression variation
├──4_pseudotime: analysis on disease trajectory
├──5_GRN: gene regulatory network inference
├──6_CCI: cell-to-cell interaction analysis
└──7_heritability: heritability analysis
```

# System requirements
Python codes require Python >= v3.8. R codes require BioC v3.17 for R >= v4.3.0.

The following version of dependencies were used when testing for compatibility.
```
anndata v0.9.1
aplot v0.2.0
Biobase v2.60.0
BiocGenerics v0.46.0
BiocParallel v1.34.2
cowplot v1.1.1
crumblr v0.99.8
dplyr v1.1.2
dreamlet v0.99.25
GenomeInfoDb v1.36.1
GenomicRanges v1.52.0
ggplot2 v3.4.3
ggtree v3.8.2
IRanges v2.34.1
limma v3.56.2
louvain v0.7.1
MatrixGenerics v1.12.3
matrixStats v1.0.0
numpy v1.24.4
pandas v1.5.0
pynndescent v0.5.6
python-igraph v0.9.10
RColorBrewer v1.1-3
S4Vectors v0.38.1
scanpy v1.9.3
scikit-learn v1.3.0
scipyv1.10.1
SingleCellExperiment v1.22.0
statsmodels v0.14.0
SummarizedExperiment v1.30.2
tidyr v1.3.0
umap v0.5.3
variancePartition v1.31.15
zellkonverter v1.10.1
zenith v1.2.0
```
For more information about the installation, demo workflow, and use cases of Dreamlet, please visit https://diseaseneurogenomics.github.io/dreamlet/ for more information.

# Dataset
Supplementary data and tables are available at [https://www.synapse.org/#!Synapse:syn52795287/wiki/624275](https://www.synapse.org/#!Synapse:syn52795287)

# Citation
Plasticity of Human Microglia and Brain Perivascular Macrophages in Aging and Alzheimer’s Disease
Lee et. al. medRxiv 2023.10.25.23297558; doi: https://doi.org/10.1101/2023.10.25.23297558

# License
MIT License

<!-- Google tag (gtag.js) -->
<script async src="https://www.googletagmanager.com/gtag/js?id=G-0CPVRBELR2"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());

  gtag('config', 'G-0CPVRBELR2');
</script>
