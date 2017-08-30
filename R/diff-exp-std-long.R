library(cudarpackage)
library(dplyr)
set.seed(24807923)
source("util/main-functions.R")

#run a pilot chain and reorder clusters
inits <- initialize_chain(sample(1e8,1), K=1.5*2^10, n.iter=200000, methodPi="Stick-breaking")
saveRDS(inits, "init-run-long-std.rds")

output <- sample_bnp_model(inits)
saveRDS(output, "long-run-std.rds")
