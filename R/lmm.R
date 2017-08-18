library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)
load(file = "diff_expr.Rdata") # loads my_data and my_data_wide; zero count genes removed

#add fixed effects
lmm_fit <- lmer(log(total+1) ~ 1 + halfdiff + (1+halfdiff|GeneID), data = my_dat)
str(lmm_fit)

library(dplyr)
wide <- Paschold2012 %>%
  mutate(genotype_replicate = paste(genotype,replicate,sep="_")) %>%
  select(GeneID, genotype_replicate, total) %>%
  tidyr::spread(genotype_replicate, total)

head(wide)

#depends on cuda_rpackage/R/data.R
source("../../../cuda_rpackage/R/data.R")
dat <- formatData(counts = my_dat_wide[,2:9], groups = rep(1:2, each=4), X = X, voom = FALSE)
est <- indEstimates(dat)

saveRDS(dat, "data/cuda-data.rds")
saveRDS(lmm_fit, "data/lmm_fit.rds")
saveRDS(est, "data/est.rds")
######## SAVED #########
