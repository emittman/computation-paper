library(cudarpackage)
library(dplyr)
set.seed(1001001001)

load("data/cuda_dat.Rdata")
load("data/ind_est.RData")


priors <- formatPriors(K=2^10, estimates = ind_est, A=3, B=3/sqrt(cuda_dat$G))

C <- list(diff_expr = matrix(c(0, 1, 0),1,3, byrow=T)) 

contr <- formatControl(n_iter = 2,
                       thin = 1,
                       warmup = 10000,
                       methodPi = "symmDirichlet",
                       idx_save = 1,
                       n_save_P = 1,
                       alpha_fixed = FALSE)

#run a pilot chain and reorder clusters
start.chain <- initFixedGrid(priors = priors, estimates = ind_est, C = C)
# init.run <- mcmc(cuda_dat, priors, contr, start.chain)
# saveRDS(init.run, "init-run.rds")
init.run <- readRDS("init-run-long.rds")
id <- order(init.run[['state']]$pi, decreasing=TRUE)
init.chain <- with(init.run[['state']], 
                   formatChain(beta = beta[,id],
                               pi = exp(pi[id]),
                               tau2 = tau2[id],
                               zeta =  start.chain$zeta,
                               alpha = alpha,
                               C = C))
contr$n_iter <- as.integer(500000)
contr$thin <- as.integer(5)
contr$warmup <- as.integer(10000)
contr$idx_save <- 0:9
contr$n_save_P <- as.integer(100)

# sb <- mcmc(cuda_dat, priors, contr, init.chain)
# saveRDS(sb, file="RDA-SD.rds")

contr$methodPi <- "stickBreaking"
sd <- mcmc(cuda_dat, priors, contr, init.chain)
saveRDS(sd, "RDA-SB-long.rds")
