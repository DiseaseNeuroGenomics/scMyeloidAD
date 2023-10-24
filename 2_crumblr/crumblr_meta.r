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

# meta
library(muscat)
library(metafor)
library(broom)
library(tidyverse)

# versions
sessionInfo()

meta_analysis = function( tabList ){
    
    # set entry names of none
    if( is.null(names(tabList)) ){
        names(tabList) = as.character(seq(length(tabList)))
    }

    # define dataset
    for( key in names(tabList) ){
        tabList[[key]]$Dataset = key
    }

    # stack datasets
    df = do.call(rbind, tabList) 

    # meta-analysis for each matching gene and assay
    # compute se from logFC and t
    # Use the fact that t = logFC / se
    df %>%
    as_tibble %>%
    group_by(assay) %>%
    do(tidy(rma( yi = logFC, sei = logFC / t, data=., method = "FE"))) %>%
    select(-term, -type) %>%
    ungroup() %>%
    mutate(FDR = p.adjust(p.value, "fdr")) %>% 
    mutate('log10FDR' = -log10(FDR))
}

plotTree = function(tree, low="grey90", mid = "red", high="darkred", xmax.scale=1.5){

    fig = ggtree(tree, branch.length = "none") + 
        geom_tiplab(color = "black", size=4, hjust=0, offset=.2) +
        theme(legend.position="top left", plot.title = element_text(hjust = 0.5))

    # get default max value of x-axis
    xmax = layer_scales(fig)$x$range$range[2]

    # increase x-axis width
    fig + xlim(0, xmax*xmax.scale) 
}

plotCoef2 = function(tab, coef, fig.tree, low="grey90", mid = "red", high="darkred", ylab){
    tab$logFC = tab$estimate
    tab$celltype = factor(tab$assay, rev(ggtree::get_taxa_name(fig.tree)))
    tab$se = tab$std.error
    fig.es = ggplot(tab, aes(celltype, logFC)) + 
        geom_hline(yintercept=0, linetype="dashed", color="grey", linewidth=1) +
        geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0) +
        # geom_point(color="dodgerblue") +
        geom_point2(aes(color=pmin(4,-log10(FDR)), size=pmin(4,-log10(FDR)))) + 
        scale_color_gradient2(name = bquote(-log[10]~FDR), limits=c(0,4), low=low, mid=mid, high=high, midpoint=-log10(0.01)) +
        scale_size_area(name = bquote(-log[10]~FDR), limits=c(0,4)) +
        geom_text2(aes(label = '+', subset=FDR < 0.05), color = "white", size=6, vjust=.3, hjust=.5) +
        theme_classic() +
        coord_flip() +
        xlab('') + 
        ylab(ylab) +
        theme(axis.text.y=element_blank(), axis.text=element_text(size = 12), axis.ticks.y=element_blank(), text = element_text(size = 20)) +
        scale_y_continuous(breaks = scales::breaks_pretty(3))
    return(fig.es)    
}

### Aging

tab.FreshMG.age = topTable(readRDS('FreshMG_cGenes.crumblr.age.sex.rds'), coef='age', number=Inf) %>% rownames_to_column('assay')
tab.PsychAD.age = topTable(readRDS('PsychAD_cGenes.crumblr.age.sex.rds'), coef='age', number=Inf) %>% rownames_to_column('assay')
res.age = meta_analysis(list(tab.FreshMG.age, tab.PsychAD.age))

tab.FreshMG.sex = topTable(readRDS('FreshMG_cGenes.crumblr.age.sex.rds'), coef='sexDiff', number=Inf) %>% rownames_to_column('assay')
tab.PsychAD.sex = topTable(readRDS('PsychAD_cGenes.crumblr.age.sex.rds'), coef='sexDiff', number=Inf) %>% rownames_to_column('assay')
res.sex = meta_analysis(list(tab.FreshMG.sex, tab.PsychAD.sex))

tab.FreshMG.age.sex = topTable(readRDS('FreshMG_cGenes.crumblr.age.sex.rds'), coef='age:sexM', number=Inf) %>% rownames_to_column('assay')
tab.PsychAD.age.sex = topTable(readRDS('PsychAD_cGenes.crumblr.age.sex.rds'), coef='age:sexM', number=Inf) %>% rownames_to_column('assay')
res.age.sex = meta_analysis(list(tab.FreshMG.age.sex, tab.PsychAD.age.sex))

### tree
pbObj_CTRL <- readRDS('FreshMG_cGenes.pbObj.rds')
hc = buildClusterTreeFromPB(pbObj_CTRL)
fig.tree = plotTree(ape::as.phylo(hc), xmax.scale=2.2) + theme(legend.position="bottom")

### effect size
fig.es1 = plotCoef2(res.age, coef='age', fig.tree, ylab='Age')
fig.es2 = plotCoef2(res.sex, coef='sexDiff', fig.tree, ylab='SexDiff')
fig.es3 = plotCoef2(res.age.sex, coef='age:sexM', fig.tree, ylab='Age:Sex')

### combine plots
options(repr.plot.width=13, repr.plot.height=5)
fig.es1 %>% insert_left(fig.tree, width=1.3) %>% insert_right(fig.es2, width=1) %>% insert_right(fig.es3, width=1)

### AD phenotypes

tab.FreshMG.dx = topTable(readRDS('FreshMG_cGenes.crumblr.dx.rds'), coef='dx_AD', number=Inf) %>% rownames_to_column('assay')
tab.PsychAD.dx = topTable(readRDS('PsychAD_cGenes.crumblr.dx.rds'), coef='dx_AD', number=Inf) %>% rownames_to_column('assay')
res.dx = meta_analysis(list(tab.FreshMG.dx, tab.PsychAD.dx))

tab.FreshMG.ce = topTable(readRDS('FreshMG_cGenes.crumblr.cerad.rds'), coef='CERAD', number=Inf) %>% rownames_to_column('assay')
tab.PsychAD.ce = topTable(readRDS('PsychAD_cGenes.crumblr.cerad.rds'), coef='CERAD', number=Inf) %>% rownames_to_column('assay')
res.ce = meta_analysis(list(tab.FreshMG.ce, tab.PsychAD.ce))

tab.FreshMG.br = topTable(readRDS('FreshMG_cGenes.crumblr.braak.rds'), coef='Braak', number=Inf) %>% rownames_to_column('assay')
tab.PsychAD.br = topTable(readRDS('PsychAD_cGenes.crumblr.braak.rds'), coef='Braak', number=Inf) %>% rownames_to_column('assay')
res.br = meta_analysis(list(tab.FreshMG.br, tab.PsychAD.br))

tab.FreshMG.de = topTable(readRDS('FreshMG_cGenes.crumblr.dementia.rds'), coef='dementia', number=Inf) %>% rownames_to_column('assay')
tab.PsychAD.de = topTable(readRDS('PsychAD_cGenes.crumblr.dementia.rds'), coef='dementia', number=Inf) %>% rownames_to_column('assay')
res.de = meta_analysis(list(tab.FreshMG.de, tab.PsychAD.de))

### tree
pbObj <- readRDS('FreshMG_cGenes.pbObj.rds')
hc = buildClusterTreeFromPB(pbObj)
fig.tree = plotTree(ape::as.phylo(hc), xmax.scale=2.2) + theme(legend.position="bottom")

### effect size
fig.es1 = plotCoef2(res.dx, coef='dx_AD', fig.tree, ylab='dx_AD')
fig.es2 = plotCoef2(res.ce, coef='CERAD', fig.tree, ylab='CERAD')
fig.es3 = plotCoef2(res.br, coef='Braak', fig.tree, ylab='Braak')
fig.es4 = plotCoef2(res.de, coef='dementia', fig.tree, ylab='Dementia')

### combine plots
options(repr.plot.width=15, repr.plot.height=5)
fig.es1 %>% insert_left(fig.tree, width=1.3) %>% insert_right(fig.es2, width=1) %>% insert_right(fig.es3, width=1) %>% insert_right(fig.es4, width=1)
