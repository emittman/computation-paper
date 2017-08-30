# s1 <- readRDS("chains/chain1_sb.rds")
# s2 <- readRDS("chains/chain2_sb.rds")
# s3 <- readRDS("chains/chain3_sb.rds")
# s4 <- readRDS("chains/chain4_sb.rds")
# class(s1) <- "myMcmcObj"
# class(s2) <- "myMcmcObj"
# class(s3) <- "myMcmcObj"
# class(s4) <- "myMcmcObj"
# out <- gelman.factors(s1,s2,s3, s4, n_iter = 50000)
# 
# library(ggplot2)
# plot.df <- data.frame(t(out[[1]]),out[[2]])
# names(plot.df) <- c("beta[1]","beta[2]","beta[3]","sigma")
# 
# library(tidyr)
# plot.df <- gather(plot.df, key=parameter, value=`hat(R)`)
# 
# ggplot(plot.df, aes(x=`hat(R)`)) + geom_histogram(bins=100) +
#   facet_wrap(~parameter, labeller=label_parsed, scales = "free")+
#   scale_x_continuous(name=expression(hat(R)))+
#   # scale_y_continuous(trans='sqrt', , breaks=c(100,signif((1:5*30)^2,2)))+
#   theme_bw()
# 
# ggsave("figures_tables/rhats-param.pdf", width=9, height=7)

#diagnostic.R in cudarpackage
gelman.factors <- function(..., n_iter){
  dims.equal <- function(...){
    identicalValue <- function(x,y) if(identical(x,y)) x else FALSE
    chains <- list(...)
    val <- Reduce(identicalValue, lapply(chains, function(ch) dim(ch$summaries$means_betas)))
    if(val[1] == FALSE) return(val) else return(TRUE)
  }
  stopifnot(all(sapply(list(...), function(ch) "myMcmcObj" %in% class(ch))))
  stopifnot(dims.equal(...))
  chains <- list(...)
  m <- length(chains)
  
  #estimate parameter variance within chain
  within.chain.beta <- lapply(chains, function(ch){
    with(ch$summaries, (meansquares_betas - means_betas^2)*n_iter/(n_iter-1))
  })
  within.chain.sigma <- lapply(chains, function(ch){
    with(ch$summaries, (meansquares_sigmas - means_sigmas^2)*n_iter/(n_iter-1))
  })
  #average within chain variance
  W.beta <- apply(plyr::laply(within.chain.beta, identity),c(2,3),mean)
  W.sigma <- apply(plyr::laply(within.chain.sigma, identity),2,mean)
  
  #estimate overall means
  ovrl.beta.mean <- apply(plyr::laply(chains, function(ch){
    ch$summaries$means_betas}),
    c(2,3),mean)
  ovrl.sigma.mean <- apply(plyr::laply(chains, function(ch){
    ch$summaries$means_sigmas}),
    2,mean)
  
  #between chain squared differences
  btween.sq.diff.beta <- lapply(chains, function(ch){
    with(ch$summaries, (means_betas - ovrl.beta.mean)^2)
  })
  btween.sq.diff.sigma <- lapply(chains, function(ch){
    with(ch$summaries, (means_sigmas - ovrl.sigma.mean)^2)
  })
  
  #between chain variances
  B.beta <- apply(plyr::laply(btween.sq.diff.beta, identity),
                  c(2,3), function(x) sum(x)*n_iter/(m-1))
  B.sigma <- apply(plyr::laply(btween.sq.diff.sigma, identity),
                   2, function(x) sum(x)*n_iter/(m-1))
  
  #estimate target variance
  targVar.beta <- (1-1/n_iter)*W.beta + 1/n_iter*B.beta
  targVar.sigma <- (1-1/n_iter)*W.sigma + 1/n_iter*B.sigma
  
  #calculate Rhat
  Rhat.beta <- sqrt(targVar.beta/W.beta)
  Rhat.sigma <- sqrt(targVar.sigma/W.sigma)
  
  list(Rhat.beta, Rhat.sigma)
  
}

combine_chains_summaries <- function(...){
  dims.equal <- function(...){
    identicalValue <- function(x,y) if(identical(x,y)) x else FALSE
    chains <- list(...)
    val <- Reduce(identicalValue, lapply(chains, function(ch) dim(ch$summaries$means_betas)))
    if(val[1] == FALSE) return(val) else return(TRUE)
  }
  stopifnot(all(sapply(list(...), function(ch) "myMcmcObj" %in% class(ch))))
  stopifnot(dims.equal(...))
  chains <- list(...)
  m <- length(chains)
  
  #average within chain variance
  beta.means <- apply(plyr::laply(chains, function(ch){ch$summaries$means_betas}),
                          c(2,3),mean)
  beta.meansquares <- apply(plyr::laply(chains, function(ch){ch$summaries$meansquares_betas}),
                                 c(2,3),mean)
  sigma.means <- apply(plyr::laply(chains, function(ch){ch$summaries$means_sigmas}),
                           2,mean)
  sigma.meansquares <- apply(plyr::laply(chains, function(ch){ch$summaries$meansquares_sigmas}),
                                  2,mean)
  probs <- apply(plyr::laply(chains, function(ch){ch$summaries$probs}, .drop = FALSE),
                      c(2,3),mean)
  list(means_betas = beta.means,
       meansquares_betas = beta.meansquares,
       means_sigmas = sigma.means,
       meansquares_sigmas = sigma.meansquares,
       prob = probs)
}
