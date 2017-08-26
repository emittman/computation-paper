initialize_chain <- function(seed, methodPi){

  cat(methodPi)  
    
  set.seed(seed)
  load("data/cuda_dat.Rdata")
  load("data/ind_est.RData")
  priors <- formatPriors(K=2^12, estimates = ind_est, A=3, B=3/sqrt(cuda_dat$G))
  
  C <- list(diff_expr = matrix(c(0, 1, 0),1,3, byrow=T)) 
  
  contr <- formatControl(n_iter = 2,
                         thin = 1,
                         warmup = 1000,
                         methodPi = "symmDirichlet",
                         idx_save = 1,
                         n_save_P = 1,
                         alpha_fixed = FALSE)
  
  #run a pilot chain and reorder clusters
  start.chain <- initFixedGrid(priors = priors, estimates = ind_est, C = C)
  init.run <- mcmc(cuda_dat, priors, contr, start.chain)
  
  id <- order(init.run[['state']]$pi, decreasing=TRUE)
  init.chain <- with(init.run[['state']],
                     formatChain(beta = beta[,id],
                                 pi = exp(pi[id]),
                                 tau2 = tau2[id],
                                 zeta =  start.chain$zeta,
                                 alpha = alpha,
                                 C = C))
  
  contr$n_iter <- as.integer(50000)
  contr$thin <- as.integer(5)
  contr$warmup <- as.integer(10000)
  contr$idx_save <- 0:9
  contr$n_save_P <- as.integer(100)
  contr$methodPi <- methodPi
  
  list(priors = priors,
       control = contr,
       seed = .Random.seed,
       init.chain = init.chain)
}

  
sample_bnp_model <- function(settings){
  load("data/cuda_dat.Rdata")
  set.seed(settings$seed)
  samples <- with(settings, mcmc(cuda_dat, priors, control, init.chain))
  samples
}
