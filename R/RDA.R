library(cudarpackage)
library(dplyr)
set.seed(1001001001)

load("data/cuda_dat.RData")
load("data/ind_est.RData")


priors <- formatPriors(K=4000, estimates = ind_est, A=3, B=3/sqrt(cuda_dat$G))

C <- list(diff_expr = matrix(c(0, 1, 0),1,3, byrow=T)) 

contr <- formatControl(n_iter = 5,
                       thin = 5,
                       warmup = 50,
                       methodPi = "stickBreaking",
                       idx_save = 1,
                       n_save_P = 1,
                       alpha_fixed = FALSE)

#run a pilot chain and reorder clusters
start.chain <- initFixedGrid(priors, ind_est)
init.run <- mcmc(cuda_dat, priors, contr, start.chain, C, ind_est)
saveRDS(init.run, "init-run.rds")
id <- order(init.run[['state']]$pi, decreasing=TRUE)
init.chain <- with(init.run[['state']], formatChain(beta[,id], exp(pi[id]), tau2[id], start.chain$zeta, alpha))
contr$n_iter <- as.integer(50000)
contr$warmup <- as.integer(1000)
contr$idx_save <- sample(cuda_dat$G, 10)
contr$n_save_P <- as.integer(100)

sb <- mcmc(cuda_dat, priors, contr, init.chain, C, ind_est)
saveRDS(sb, file="RDA-SB.rds")

contr$methodPi <- "symmDirichlet"
sd <- mcmc(cuda_dat, priors, contr, C=C, estimates=ind_est)
saveRDS(s, "RDA-SD.rds")
