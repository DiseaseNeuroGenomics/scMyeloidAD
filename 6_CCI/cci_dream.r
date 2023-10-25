#' Perform Dream linear mixed effects regression
#' on LIANA cell cell interaction scores 

library(tidyverse)
library(variancePartition)

# load metadata
meta = readRDS(metadata_path)

# load merged liana scores  
df_cci = readRDS(merged_path)

# regression formula
form = ~ dx + covariates

# make contrast
L = makeContrastsDream(form, meta, contrasts=c(AD_CTRL = "dxAD - dxCTRL"))

df_cci = -log10(df_cci) # log transform and invert magnitude_rank scores

variances <- apply(df_cci, 1, var) # use regular var
df_cci = df_cci[variances > 0,]

n_cores = 8 # number of cores for multithreading

param = SnowParam(n_cores, "SOCK", progressbar=TRUE)

fit = dream(df_cci,
            form, meta, L=L,
            BPPARAM=param
           )

# ebayes
fit = eBayes(fit)

fit_table = topTable(fit,
                     coef="AD_CTRL",
                     number = Inf,
                     sort.by = "logFC")

fit_table %>%
arrange(desc(logFC))

write.csv(fit_table, '/path/to/output.csv', row.names=TRUE)
