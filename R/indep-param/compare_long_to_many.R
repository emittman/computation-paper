source("R/util/rhat.R")

ch1 <- readRDS("chains/chain1_sb.rds")
ch2 <- readRDS("chains/chain2_sb.rds")
ch3 <- readRDS("chains/chain3_sb.rds")
ch4 <- readRDS("chains/chain4_sb.rds")
class(ch1) <- "myMcmcObj"
class(ch2) <- "myMcmcObj"
class(ch3) <- "myMcmcObj"
class(ch4) <- "myMcmcObj"

to.compare <- combine_chains_summaries(ch1,ch2,ch3,ch4)

long.chain <- readRDS("R/RDA-SB-long.rds")

plot(to.compare$means_betas[1,], long.chain$summaries$means_betas[1,])
hist(abs(to.compare$means_betas[1,] - long.chain$summaries$means_betas[1,]), 100)
hist(abs(to.compare$means_betas[2,] - long.chain$summaries$means_betas[2,]), 100)
hist(abs(to.compare$means_betas[3,] - long.chain$summaries$means_betas[3,]), 100)

library(coda)

alpha <- mcmc(long.chain$samples$alpha)
traceplot(alpha)
effectiveSize(alpha)

num.occ <- mcmc(long.chain$samples$num_occupied[1:1000*100])
traceplot(num.occ)

max.id <- mcmc(long.chain$samples$max_id[1:1000*100])
traceplot(max.id)
