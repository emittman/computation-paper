library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)
load(file = "data/my_dat.RData")
head(my_dat)

#Fit linear mixed-effects model
lmm_fit <- lmer(total ~ 1 + halfdiff + flow_cell + (1+flow_cell+halfdiff|GeneID), data = my_dat)

save(lmm_fit, file="data/lmm_fit.RData")

