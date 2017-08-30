library(cudarpackage)
library(dplyr)
set.seed(24807923)
source("../util/main-functions.R")
load("/data/cuda-dat-ind.Rdata")

#run a pilot chain and reorder clusters
inits <- initialize_chain(sample(1e8,1), K=2^10, n.iter=100, methodPi="Stick-breaking")

inits1 <- inits
inits1$control$idx_save <- as.integer(0:(cuda_dat$G-1))

inits2 <- inits
inits2$control$idx_save <- as.integer(0)

output1 <- sample_bnp_model(inits1)

output2 <- sample_bnp_model(inits2)

saveRDS(output1, "bottleneck-100iter.rds")
saveRDS(output2, "no-bottleneck-100iter.rds")
