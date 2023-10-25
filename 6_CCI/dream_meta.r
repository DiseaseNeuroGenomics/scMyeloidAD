# Dream meta analysis for FreshMG/PsychAD
# Differentially expressed cell-cell interactions

library(metafor)
library(tidyverse)
library(broom)

fmg_table = read.csv('/path/to/FreshMG_toptable.csv')
pad_table = read.csv('/path/to/PsychAD_toptable.csv')

meta_analysis = function( tabList ){
    options(na.action = "na.pass")
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

    # make CCI column from index 
    # and drop prepended dataset number
    df['CCI'] = row.names(df) %>% str_sub(start=3)
    
    # meta-analysis for each matching gene and assay
    # compute se from logFC and t
    # Use the fact that t = logFC / se
     df %>% 
        as_tibble %>%
        group_by(CCI) %>% # group by CCI column
        do(tidy(rma( yi = logFC, sei = logFC / t, data=., method = "FE"))) %>%
        select(-term, -type)
}

# run analysis
tabList = list(fmg_table, pad_table)

res = meta_analysis( tabList )

res$adj.P.Val = p.adjust(res$p.value, "fdr")

write.csv(res, '/path/to/output.csv', row.names=TRUE)