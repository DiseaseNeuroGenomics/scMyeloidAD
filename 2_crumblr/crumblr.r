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

proc_crumblr = function(pbObj, form, coef, save=NA){

    ### analysis with dream()
    cobj = crumblr(cellCounts(pbObj))
    fit = dream(cobj, form, colData(pbObj))
    fit = eBayes(fit)
    
    if(!is.na(save)){
        saveRDS(fit, save)
    }
    
    ### Multivariate test along hierarchy
    hc = buildClusterTreeFromPB(pbObj)
    res = treeTest(fit, cobj, hc, coef=coef)
    
    ### tree
    fig.tree = plot_tree(res, xmax.scale=1.4) + theme(legend.position="bottom")

    ### effect size
    fig.es = plot_es(fit, coef, fig.tree)
    
    ### Variance partitioning analysis
    fig.pct = plot_vp_pct(pbObj, form)
    
    ### combine plots
    options(repr.plot.width=13, repr.plot.height=5)
    fig.es %>% insert_left(fig.tree, width=4) %>% insert_right(fig.pct, width=1)
}

plot_vp = function(pbObj, form){
    
    ### Variance partitioning analysis
    cobj = crumblr(cellCounts(pbObj))
    vp.c = fitExtractVarPartModel(cobj, form, colData(pbObj))

    fig.vp = plotVarPart(sortCols(vp.c), label.angle=60, ncol=4) + theme(aspect.ratio=1)
    return(fig.vp)
}

plot_vp_pct = function(pbObj, form){
    
    ### Variance partitioning analysis
    cobj = crumblr(cellCounts(pbObj))
    vp.c = fitExtractVarPartModel(cobj, form, colData(pbObj))

    fig.pct = plotPercentBars( sortCols(vp.c) ) + 
              theme(legend.position="right", text = element_text(size = 12), axis.text.y=element_blank(), axis.ticks.y=element_blank())

    return(fig.pct)
}

plot_tree = function(tree, low="grey90", mid = "red", high="darkred", xmax.scale=1.5){

    # PASS R check
    isTip = label = node = FDR = NULL

    fig.tree = ggtree(tree, branch.length = "none") + 
               geom_tiplab(color = "black", size=4, hjust=0, offset=.2) +
               geom_point2(aes(label = node, color=pmin(4,-log10(FDR)), size=pmin(4,-log10(FDR)))) + 
               scale_color_gradient2(name = bquote(-log[10]~FDR), limits=c(0,4), low=low, mid=mid, high=high, midpoint=-log10(0.01)) +
               scale_size_area(name = bquote(-log[10]~FDR), limits=c(0,4)) +
               geom_text2(aes(label = '+', subset=FDR < 0.05), color = "white", size=6, vjust=.3, hjust=.5) +
               theme(legend.position="top left", plot.title = element_text(hjust = 0.5))

    # get default max value of x-axis
    xmax = layer_scales(fig.tree)$x$range$range[2]

    # increase x-axis width
    fig.tree = fig.tree + xlim(0, xmax*xmax.scale)
    
    return(fig.tree)
}

plot_tree_simple = function(tree, low="grey90", mid = "red", high="darkred", xmax.scale=1.5){

    fig.tree = ggtree(tree, branch.length = "none") + 
               geom_tiplab(color = "black", size=4, hjust=0, offset=.2) +
               theme(legend.position="top left", plot.title = element_text(hjust = 0.5))

    # get default max value of x-axis
    xmax = layer_scales(fig.tree)$x$range$range[2]

    # increase x-axis width
    fig.tree = fig.tree + xlim(0, xmax*xmax.scale) 
    
    return(fig.tree)
}

plot_es = function(fit, coef, tree){
    tab = topTable(fit, coef=coef, number=Inf)
    tab$celltype = factor(rownames(tab), rev(get_taxa_name(tree)))
    tab$se = with(tab, logFC/t)
    fig.es = ggplot(tab, aes(celltype, logFC)) + 
      geom_hline(yintercept=0, linetype="dashed", color="grey", linewidth=1) +
      geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0) +
      geom_point(color="dodgerblue") +
      theme_classic() +
      coord_flip() +
      xlab('') + 
      ylab("Effect size") +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), text = element_text(size = 12))
    return(fig.es)
}

plot_es_coef = function(fit, coef, tree, low="grey90", mid = "red", high="darkred"){
    ### effect size
    tab = topTable(fit, coef=coef, number=Inf)
    tab$celltype = factor(rownames(tab), rev(get_taxa_name(tree)))
    tab$se = with(tab, logFC/t)
    fig.es = ggplot(tab, aes(celltype, logFC)) + 
        geom_hline(yintercept=0, linetype="dashed", color="grey", linewidth=1) +
        geom_errorbar(aes(ymin = logFC - 1.96*se, ymax = logFC + 1.96*se), width=0) +
        # geom_point(color="dodgerblue") +
        geom_point2(aes(color=pmin(4,-log10(adj.P.Val)), size=pmin(4,-log10(adj.P.Val)))) + 
        scale_color_gradient2(name = bquote(-log[10]~adj.P.Val), limits=c(0,4), low=low, mid=mid, high=high, midpoint=-log10(0.01)) +
        scale_size_area(name = bquote(-log[10]~adj.P.Val), limits=c(0,4)) +
        geom_text2(aes(label = '+', subset=adj.P.Val < 0.05), color = "white", size=6, vjust=.3, hjust=.5) +
        theme_classic() +
        coord_flip() +
        xlab('') + 
        ylab(coef) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), text = element_text(size = 12))
    return(fig.es)    
}

prefix = 'example_run'
form = ~ log(n_counts) + (1|batch) + scale(age) + (1|sex) + (1|ancestry) + PMI + Braak + 1

### loading pbObj
pbObj <- readRDS(paste0(prefix,'.pbObj.rds'))

### Variance partitioning analysis
options(repr.plot.width=6, repr.plot.height=6)
plot_vp(pbObj, form = form)

### crumblr
proc_crumblr(pbObj, form = form, coef = 'Braak', save=paste0(prefix,'.crumblr.braak.rds'))
