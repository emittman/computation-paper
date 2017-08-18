#expected no. of clusters given alpha and G
exp_Kocc <- function(alpha, G){
  sum(alpha/(alpha + 1:G - 1))
}

#expected no. clusters (expressed as a power of G)
# given G and alpha (expressed as a power of G)
power.of.G <- Vectorize(
  function(G, pow){
  log(exp_Kocc(G^pow, G))/log(G)
  }, vectorize.args = "G"
)

x <- exp(seq(log(5000), log(40000), length.out=20))
plot(x,power.of.G(x, .3))

#interval for no. of clusters given gamma prior for alpha in terms
# of expected value as a power of G and shape parameter
exp_Kocc_gam <- function(pow, shape, G, q=c(.01,.5,.99)){
  qg <- qgamma(q, shape, shape/G^pow)
  sapply(qg, function(x) exp_Kocc(x, G))
}

gam_implied_alpha <- function(G, pow, shape, q=c(.01,.5,.99)){
  log(qgamma(q, shape, shape/G^pow))/log(G)
}

#interval for no. of clusters (expressed as a power of G)
# given gamma prior for alpha in terms of expected value as a power of G
# and a scale parameter
power.of.G_gam <- function(G, pow, shape, q=c(.01,.5,.99)){
  k <- exp_Kocc_gam(pow, shape, G, q)
  log(k)/log(G)
}

exp_Kocc_gam(.3, 3, 40000, c(.001,.999))
power.of.G_gam(40000, .3, 3)

library(plyr)
library(ggplot2)
x <- exp(seq(log(5000), log(40000), length.out=5))
df <- ldply(x, function(i){
  sdf <- data.frame(t(power.of.G_gam(i, .5, 3)))
  names(sdf) <- c("lower01", "median", "upper99")
  sdf
})

p1 <- ggplot(df, aes(x=x, y=median, ymin=lower01, ymax=upper99))+
  geom_pointrange() +
  scale_x_continuous(trans="log", breaks=signif(x,2))+
  ylab("Expected no. occ. as power of G") +
  xlab("")+
  ggtitle(expression(paste(alpha %~% Gamma(3,3/G^.5)))) +
  theme_classic(base_size = 12)
df2 <- ldply(x, function(i){
  sdf <- data.frame(t(exp_Kocc_gam(.5, 3, i)))
  names(sdf) <- c("lower01", "median", "upper99")
  sdf
})
p2 <- ggplot(df2, aes(x=x, y=median, ymin=lower01, ymax=upper99))+
  geom_pointrange() +
  scale_x_continuous(trans="log", breaks=signif(x,2))+
  ylab("Expected no. occ.") + xlab("G")+
  theme_classic(base_size=12)
library(cowplot)
pdf("figures_tables/alphaprior.pdf")
plot_grid(p1, p2, ncol=1)
dev.off()
