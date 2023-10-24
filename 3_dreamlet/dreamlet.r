# dreamlet universe
library(dreamlet)
library(crumblr)

# data IO
library(SingleCellExperiment) 
library(zellkonverter)
library(tidyr)

# plotting
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(ggtree)
library(aplot)

# versions
sessionInfo()

prefix = 'example_run'
form = ~ log(n_counts) + (1|batch) + scale(age) + (1|sex) + (1|ancestry) + PMI + Braak + 1

### loading pbObj
pbObj <- readRDS(paste0(prefix,'.pbObj.rds'))

### Normalize and apply voom/voomWithDreamWeights
res.proc = processAssays(pbObj, form, min.cells=5, min.count=5, min.samples=4, min.prop=0.2)
saveRDS(res.proc, paste0(prefix,'.res.proc.rds'))

# variance partition
vp.lst = fitVarPart(res.proc, form)
saveRDS(vp.lst, paste0(prefix,'.vp.lst.rds'))

# number of genes in VP
length(unique(vp.lst$gene))

# aggregate
vp.agg <- aggregate(. ~ gene, vp.lst[,-1], mean)
rownames(vp.agg) <- vp.agg$gene
vp.agg <- vp.agg[,-1]

options(repr.plot.width=5, repr.plot.height=5)
plotVarPart(sortCols(vp.agg), label.angle=45, ncol=1) + theme(aspect.ratio=1)

### Differential expression analysis within each assay
res.dl.braak = dreamlet(res.proc, form)
saveRDS(res.dl.braak, paste0(prefix,'.res.dl.braak.rds'))
