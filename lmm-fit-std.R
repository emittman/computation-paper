library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)
load(file = "data/my_dat_std.RData")
head(my_dat)

#Fit linear mixed-effects model
lmm_fit <- lmer(total ~ 1 + genotype + flow_cell + (1 + genotype + flow_cell|GeneID), data = my_dat,
                control = lmerControl(optCtrl = list(maxfun=1000)))

max(abs(with(lmm_fit@optinfo$derivs,solve(Hessian,gradient))))

save(lmm_fit, file="data/lmm_fit_std.RData")
