load("data/my_dat_wide_std.RData")
bnp.fit <- readRDS("long-run-std.rds")
bnp <- with(bnp.fit$summaries, data.frame(t(rbind(means_betas,means_sigmas))))
names(bnp) <- c("intercept","Mo17","flow_cell","sigma")
bnp$GeneID <- my_dat_wide$GeneID
bnp_long <- tidyr::gather(bnp, key=genotype, value=est, -GeneID)
bnp_long$type <- "BNP"
library(dplyr)
bnp.top <- arrange(bnp, abs(Mo17)) %>% select(GeneID) %>% head(35000)

load("data/ind-est-std.RData")
ols <- data.frame(t(ind_est$beta))
names(ols) <- c("intercept","Mo17","flow_cell")
ols$GeneID <- my_dat_wide$GeneID
ols_long <- tidyr::gather(ols, key=genotype, value=est, -GeneID)
ols_long$type <- "ols"



load("data/lmm_fit_std.RData")
lmm <- data.frame(coef(lmm_fit)[[1]])
names(lmm) <- c("intercept","Mo17","flow_cell")
lmm$GeneID <- my_dat_wide$GeneID
lmm.top <- arrange(lmm, -abs(Mo17)) %>% select(GeneID) %>% head(1000)

which(lmm.top[[1]] %in% bnp.top[[1]])
#[1] 86, 953
#id <- with(ols, which(Mo17< -2 & Mo17 > -2.1& intercept<5 & intercept>4.9))
# id <- with(ols, which(abs(Mo17)>1 & abs(Mo17)<1.1 & sigma>.25 & sigma<.28))
id <- with(ols, which(abs(flow_cell)>3))
id2 <- with(ols, which(flow_cell>.4 & flow_cell<.42))


library(ggplot2)
bnp$method="bnp"
ols$method="ols"
lmm$method="lmm"
df <- rbind(bnp[,-4],ols,lmm)



# id <- which(lmm$GeneID == lmm.top$GeneID[720])
bnp[id,]
lmm[id,]

my_dat_wide[id,]
with(ols, plot(Mo17, flow_cell, pch="."))#, xlim=c(-2.2, -1.8), ylim=c(-.2,.5)))
with(ols, points(Mo17[id], flow_cell[id], pch=4, col=4, cex=1.2))
with(lmm, points(Mo17[id], flow_cell[id], pch=2,col=2, cex=1.2))
with(bnp, points(Mo17[id], flow_cell[id], pch=3,col=3, cex=1.2))
abline(h=0)
abline(v=0)
legend(4, -1.5, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)


with(ols, plot(intercept, Mo17, pch=".", ylim=c(-.1,1.5)))
with(ols, points(intercept[id], Mo17[id], pch=4, col=4))
with(lmm, points(intercept[id], Mo17[id], pch=2,col=2))
with(bnp, points(intercept[id], Mo17[id], pch=3,col=3))
legend(8, 1.5, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)

with(ols, plot(intercept, flow_cell, pch="."))
with(ols, points(intercept[id], flow_cell[id], pch=4, col=4))
with(lmm, points(intercept[id], flow_cell[id], pch=2,col=2))
with(bnp, points(intercept[id], flow_cell[id], pch=3,col=3))
legend(8, -2, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)

ols$sigma <- sqrt(ind_est$sigma2-.1)
with(ols, plot(intercept, sigma, pch="."))
with(ols, points(intercept[id], sigma[id], pch=4, col=4))
with(bnp, points(intercept[id], sigma[id], pch=3,col=3))
legend(6, 1.2, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)

with(ols, plot(Mo17, sigma, pch="."))
with(ols, points(Mo17[id], sigma[id], pch=4, col=4))
with(bnp, points(Mo17[id], sigma[id], pch=3,col=3))
legend(-6, 1.2, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)

#id2
with(ols, plot(Mo17, flow_cell, pch="."))#, xlim=c(-2.2, -1.8), ylim=c(-.2,.5)))
with(ols, points(Mo17[id2], flow_cell[id2], pch=4, col=4, cex=1.2))
with(lmm, points(Mo17[id2], flow_cell[id2], pch=2,col=2, cex=1.2))
with(bnp, points(Mo17[id2], flow_cell[id2], pch=3,col=3, cex=1.2))
abline(h=0)
legend(4, -1.5, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)


with(ols, plot(intercept, Mo17, pch=".", ylim=c(-.1,1.5)))
with(ols, points(intercept[id2], Mo17[id2], pch=4, col=4))
with(lmm, points(intercept[id2], Mo17[id2], pch=2,col=2))
with(bnp, points(intercept[id2], Mo17[id2], pch=3,col=3))
legend(8, 1.5, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)

with(ols, plot(intercept, flow_cell, pch="."))
with(ols, points(intercept[id2], flow_cell[id2], pch=4, col=4))
with(lmm, points(intercept[id2], flow_cell[id2], pch=2,col=2))
with(bnp, points(intercept[id2], flow_cell[id2], pch=3,col=3))
legend(8, -2, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)

ols$sigma <- sqrt(ind_est$sigma2-.1)
with(ols, plot(intercept, sigma, pch="."))
with(ols, points(intercept[id2], sigma[id2], pch=4, col=4))
with(bnp, points(intercept[id2], sigma[id2], pch=3,col=3))
legend(6, 1.2, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)

with(ols, plot(Mo17, sigma, pch="."))
with(ols, points(Mo17[id2], sigma[id2], pch=4, col=4))
with(bnp, points(Mo17[id2], sigma[id2], pch=3,col=3))
legend(-6, 1.2, legend=c("ols", "lmm", "bnp"),
       col=c("blue", "red", "green"), pch=c(4,2,3), cex=0.8)
