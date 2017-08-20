library(cudarpackage)
library(dplyr)
set.seed(1001001001)

load("diff_expr.Rdata")
dat <- readRDS("data/cuda-data.rds")
est <- readRDS("data/est.rds")

X <- matrix(c(1,  -1,
              1,   1),byrow=T,2,2)


priors <- formatPriors(K=4000, estimates = est, A=3, B=3/sqrt(dat$G))

C <- list(diff_expr = matrix(c(0, 1),1,2, byrow=T)) 

contr <- formatControl(n_iter = 5,
                       thin = 5,
                       warmup = 1000,
                       methodPi = "stickBreaking",
                       idx_save = 1,
                       n_save_P = 1,
                       alpha_fixed = FALSE)

#run a pilot chain and reorder clusters
start.chain <- initFixedGrid(priors, est)
init.run <- mcmc(dat, priors, contr, start.chain, C, est)
id <- order(init.run[['state']]$pi, decreasing=TRUE)
init.chain <- with(init.run[['state']], formatChain(beta[,id], exp(pi[id]), tau2[id], zeta, alpha))
contr$n_iter <- 100000
contr$warmup <- 10000
contr$idx_save <- sample(dat$G, 10)
contr$n_save_P <- 100

sb <- mcmc(dat, priors, contr, init.chain, C, est)
saveRDS(sb, file="RDA-SB.rds")

contr$methodPi <- "symmDirichlet"
sd <- mcmc(dat, priors, contr, C=C, estimates=est)
saveRDS(s, "RDA-SB.rds")
