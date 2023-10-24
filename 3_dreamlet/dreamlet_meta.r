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
    group_by(ID, assay) %>%
    do(tidy(rma( yi = logFC, sei = logFC / t, data=., method = "FE"))) %>%
    select(-term, -type) %>%
    ungroup() %>%
    group_by(assay) %>%
    mutate(FDR = p.adjust(p.value, "fdr")) %>%
    mutate('log10FDR' = -log10(FDR))
}

prefix = 'example_run'

# run meta analysis
tab.FreshMG.ag = topTable(readRDS('FreshMG_cGenes.res.dl.age.rds'), coef='scale(age)', number=Inf)
tab.PsychAD.ag = topTable(readRDS('PsychAD_cGenes.res.dl.age.rds'), coef='scale(age)', number=Inf)
res.ag = meta_analysis(list(tab.FreshMG.ag, tab.PsychAD.ag))
saveRDS(res.ag, paste0(prefix,'.res.ag.rds'))

tab.FreshMG.dx = topTable(readRDS('FreshMG_cGenes.res.dl.dx.rds'), coef='dx_AD', number=Inf)
tab.PsychAD.dx = topTable(readRDS('PsychAD_cGenes.res.dl.dx.rds'), coef='dx_AD', number=Inf)
res.dx = meta_analysis(list(tab.FreshMG.dx, tab.PsychAD.dx))
saveRDS(res.dx, paste0(prefix,'.res.dx.rds'))

tab.FreshMG.ce = topTable(readRDS('FreshMG_cGenes.res.dl.cerad.rds'), coef='CERAD', number=Inf)
tab.PsychAD.ce = topTable(readRDS('PsychAD_cGenes.res.dl.cerad.rds'), coef='CERAD', number=Inf)
res.ce = meta_analysis(list(tab.FreshMG.ce, tab.PsychAD.ce))
saveRDS(res.ce, paste0(prefix,'.res.ce.rds'))

tab.FreshMG.br = topTable(readRDS('FreshMG_cGenes.res.dl.braak.rds'), coef='Braak', number=Inf)
tab.PsychAD.br = topTable(readRDS('PsychAD_cGenes.res.dl.braak.rds'), coef='Braak', number=Inf)
res.br = meta_analysis(list(tab.FreshMG.br, tab.PsychAD.br))
saveRDS(res.br, paste0(prefix,'.res.br.rds'))

tab.FreshMG.de = topTable(readRDS('FreshMG_cGenes.res.dl.dementia.rds'), coef='dementia', number=Inf)
tab.PsychAD.de = topTable(readRDS('PsychAD_cGenes.res.dl.dementia.rds'), coef='dementia', number=Inf)
res.de = meta_analysis(list(tab.FreshMG.de, tab.PsychAD.de))
saveRDS(res.de, paste0(prefix,'.res.de.rds'))

