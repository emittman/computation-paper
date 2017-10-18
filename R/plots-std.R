library(dplyr)
library(lme4)
library(ggplot2)
library(GGally)

load("data/lmm_fit_std.RData")
load("data/ind-est-std.RData")
load("data/my_dat_wide_std.RData")

ols <- data.frame(t(ind_est$beta))
names(ols) <- c("intercept","Mo17","flow_cell")
ols$GeneID <- my_dat_wide$GeneID
ols_long <- tidyr::gather(ols, key=genotype, value=est, -GeneID)
ols_long$type <- "independent"

lmm <- coef(lmm_fit)[[1]]
names(lmm) <- c("intercept","Mo17","flow_cell")
lmm$GeneID <- my_dat_wide$GeneID
lmm_long <- tidyr::gather(lmm, key=genotype, value=est, -GeneID)
lmm_long$type <- "hierarchical"

bnp.fit <- readRDS("long-run-std.rds")
bnp <- with(bnp.fit$summaries, data.frame(t(rbind(means_betas,means_sigmas))))
names(bnp) <- c("intercept","Mo17","flow_cell","sigma")
bnp$GeneID <- my_dat_wide$GeneID
bnp_long <- tidyr::gather(bnp, key=genotype, value=est, -GeneID)
bnp_long$type <- "BNP"

source("util/rhat.R")
bnp.chain1 <- readRDS("chains/chain1_sb.rds"); class(bnp.chain1) <- "myMcmcObj"
bnp.chain2 <- readRDS("chains/chain2_sb.rds"); class(bnp.chain2) <- "myMcmcObj"
bnp.chain3 <- readRDS("chains/chain3_sb.rds"); class(bnp.chain3) <- "myMcmcObj"
bnp.chain4 <- readRDS("chains/chain4_sb.rds"); class(bnp.chain4) <- "myMcmcObj"
bnp.fit2 <- combine_chains_summaries(bnp.chain1,bnp.chain2,bnp.chain3)
bnp2 <- with(bnp.fit2, data.frame(t(rbind(means_betas, log(means_sigmas)))))
names(bnp2) <- c("intercept","Mo17","flow_cell","log_sigma")
bnp2$GeneID <- my_dat_wide$GeneID
bnp_long <- tidyr::gather(bnp2, key=genotype, value=est, -GeneID)
bnp_long$type <- "BNP"

id <- sample(my_dat_wide$GeneID, size=40)

rhat_df <- gelman.factors(bnp.chain1,bnp.chain2,bnp.chain3,bnp.chain4, n_iter = 1000)
rhat_df <- data.frame(cbind(t(rhat_df[[1]]),rhat_df[[2]]))
names(rhat_df) <- c("intercept","Mo17","flow cell", "sigma")
gather(rhat_df, key=par, value=value) %>%
  ggplot(aes(value)) + geom_histogram(bins=100) + facet_wrap(~par, scales = "free") + theme_bw()
ggsave("../figures_tables/rhats-param.pdf", width=5, height=5)


plot_df <- rbind(filter(ols_long),# GeneID %in% id),
                 filter(lmm_long),
                 filter(bnp_long,genotype!='log_sigma'))#, GeneID %in% id))


library(tidyr)
my_theme <- theme_bw(base_size=12)
no_panel <- theme(panel.grid = element_blank())
no_axes <-  theme(axis.text = element_blank(),
                  axis.ticks = element_blank())

#figure 1
my_hex <- function(data, mapping, ...){
  ggplot(data, mapping, ...) + geom_hex(bins=20) +
    scale_fill_continuous(low="white",high="darkblue",trans="log",
                          breaks=c(1,10,100,1000,10000))+
    my_theme + no_panel
}
my_dens <- function(data, mapping, ...){
  ggplot(data, mapping, ...) + geom_density(fill="darkblue")+
    my_theme + no_panel
}
mutate(ols, sigma=sqrt(ind_est$sigma2)) %>%
  select(intercept, Mo17, flow_cell, sigma) %>%
  ggpairs(columnLabels = c(paste(expression(beta[1])),
                           paste(expression(beta[2])),
                           paste(expression(beta[3])),
                           paste(expression(sigma))),
          lower=list(continuous=my_hex),
          upper=list(continuous=my_hex),
          diag=list(continuous=my_dens), labeller = label_parsed)
ggsave("../figures_tables/pairs-ind-est-std.pdf", width=9, height=8.25)

#####
# Compare scalar components by method


plot_df2 <- spread(plot_df, key=type, value=est)
# plot_df2$genotype <- ifelse(plot_df2$genotype=="intercept", paste(expression(hat(beta)[1])),
                                     # ifelse(plot_df2$genotype=="Mo17", paste(expression(hat(beta)[2])),
                                            # paste(expression(hat(beta)[3]))))
plot_df2$genotype <- factor(plot_df2$genotype, levels=c("intercept","Mo17","flow_cell"),
                            labels = c("intercept", "Mo17","flow cell"))
p1 <- ggplot(plot_df2, aes(x=ols, y=lmm)) + geom_hex(bins=30) + 
  geom_abline(slope=1,alpha=.2) + 
  facet_wrap(~genotype, scales="free")+#, labeller=label_parsed) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))+
  my_theme
#+
  # xlab("OLS") + ylab("LMM")

p1.1 <- ggplot(plot_df2, aes(x=ols, y=BNP)) + geom_hex(bins=30) + 
  geom_abline(slope=1,alpha=.2) + 
  facet_wrap(~genotype, scales="free")+#, labeller=label_parsed) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))+
  my_theme#+
  # xlab("OLS") + ylab("BNP")

p1.2 <- ggplot(plot_df2, aes(x=lmm, y=BNP)) + geom_hex(bins=30) + 
  geom_abline(slope=1,alpha=.2) + 
  facet_wrap(~genotype, scales="free")+#, labeller=label_parsed) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))+
  xlab("hierarchical") + ylab("BNP")

# leg <- cowplot::get_legend(p1)
# tm.noleg <- theme(legend.position = "none")
# plt.grd1 <- cowplot::plot_grid(cowplot::plot_grid(p1+tm.noleg, p1.1+tm.noleg, ncol=1),leg,rel_widths = c(1,.15))

ggmatrix(list(p1,p1.1), 2,1,yAxisLabels = c("hierarch. vs. indep.", "BNP vs. indep."))

ggsave("../figures_tables/method-compare-std.pdf", width=10, height=6)

#################
# Pairwise joint estimates
plot_df3 <- spread(plot_df, key=genotype, value=est)


# p2 <- ggplot(plot_df3, aes(x=intercept, y=halfdiff)) + geom_hex() +
#   facet_wrap(~type, labeller = as_labeller(facet_names))+
#   scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
#   theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
#                                strip.text.x = element_text(size=12))



# plot_df_ws <- cbind(ols, log_sigma = .5*log(ind_est$sigma2)) %>%
#   mutate(GeneID = my_dat_wide$GeneID)

library(GGally)

### Need to set common scales manually
my_x_axis <-  geom_hline(yintercept=0)
my_y_axis <- geom_vline(xintercept=0)
q11 <- filter(plot_df3, type=="independent") %>%
  ggplot(aes(x=intercept, y=Mo17)) + geom_hex(bins=20) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col1xlim <- ggplot_build(q11)$layout$panel_ranges[[1]]$x.range
col1ylim <- ggplot_build(q11)$layout$panel_ranges[[1]]$y.range


q12 <- filter(plot_df3, type=="independent") %>%
  ggplot(aes(x=intercept, y=flow_cell)) + geom_hex(bins=20) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col2xlim <- ggplot_build(q12)$layout$panel_ranges[[1]]$x.range
col2ylim <- ggplot_build(q12)$layout$panel_ranges[[1]]$y.range


q13 <- filter(plot_df3, type=="independent") %>%
  ggplot(aes(x=Mo17, y=flow_cell)) + geom_hex(bins=20) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col3xlim <- ggplot_build(q13)$layout$panel_ranges[[1]]$x.range
col3ylim <- ggplot_build(q13)$layout$panel_ranges[[1]]$y.range
# q13 <- q13 + ylim(col3ylim) + xlim(col3xlim)

q21 <- filter(plot_df3, type=="hierarchical") %>%
  ggplot(aes(x=intercept, y=Mo17)) + geom_hex(bins=20) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme + xlim(col1xlim) + ylim(col1ylim)+theme(axis.text = element_blank(),
                                                   axis.ticks = element_blank())

q22 <- filter(plot_df3, type=="hierarchical") %>%
  ggplot(aes(x=intercept, y=flow_cell)) + geom_hex(bins=20) + 
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme + xlim(col2xlim) + ylim(col2ylim)+theme(axis.text = element_blank(),
                                                   axis.ticks = element_blank())
# #geom_violin(color="darkblue", alpha=0)+
# geom_vline(xintercept=log(summary(lmm_fit)$sigma), linetype=2, color="darkblue")+
# my_theme + coord_flip() + ylim(col2xlim) + xlim(col2ylim)

q23 <- filter(plot_df3, type=="hierarchical") %>%
  ggplot(aes(x=Mo17, y=flow_cell)) + geom_hex(bins=20) + 
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme + xlim(col3xlim) + ylim(col3ylim)

q31 <- bnp2 %>%
  ggplot(aes(x=intercept, y=Mo17)) + geom_hex(bins=20) +
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme + xlim(col1xlim) + ylim(col1ylim)

q32 <- bnp2 %>%
  ggplot(aes(x=intercept, y=flow_cell)) + geom_hex(bins=20) +
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme + xlim(col2xlim) + ylim(col2ylim) + theme(axis.text = element_blank(),
                                                      axis.ticks = element_blank())

q33 <- bnp2 %>%
  ggplot(aes(x=Mo17, y=flow_cell)) + geom_hex(bins=20) +
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme + xlim(col3xlim) + ylim(col3ylim)

qqq1<- ggmatrix(list(q11,q21,q31,
                     q12,q22,q32,
                     q13,q23,q33),3,3,xAxisLabels = c("independent","hierarchical","BNP"),
                yAxisLabels = c("Mo17 vs. intercept",
                                "flow cell vs. intercept",
                                "flow cell vs. Mo17"),
                # yAxisLabels = c(paste(expression(paste(hat(beta)[2],"  vs.  ",hat(beta)[1]))),
                                # paste(expression(paste(hat(beta)[3],"  vs.  ",hat(beta)[1]))),
                                # paste(expression(paste(hat(beta)[3],"  vs.  ",hat(beta)[2])))),
                # labeller = label_parsed, legend=c(1,1),
                showStrips = TRUE) + 
  theme(strip.text.y = element_text(size=12),
        strip.text.x = element_text(size=12))

ggsave("../figures_tables/pairs-3-methods-std.pdf", width=10, height=9)


#mean-variance trend log(sigma) vs beta[1]
mv.trend.df <- data.frame(rbind(cbind(ols[,1:3], sigma = sqrt(ind_est$sigma2), type="independent"),
                                cbind(bnp2[,1:3], sigma = exp(bnp2$log_sigma), type="BNP"))
)
mv11 <- filter(mv.trend.df, type=="independent") %>%
  ggplot(aes(intercept, sigma)) + geom_hex(bins=20) + facet_wrap(~type)+
    scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
    my_theme

row1x <- ggplot_build(mv11)$layout$panel_ranges[[1]]$x.range
row1y <- ggplot_build(mv11)$layout$panel_ranges[[1]]$y.range

mv12 <- filter(mv.trend.df, type=="BNP") %>%
  ggplot(aes(intercept, sigma)) + geom_hex(bins=20) + facet_wrap(~type)+
    scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
    my_theme + xlim(row1x) + ylim(row1y)


mv21 <- filter(mv.trend.df, type=="independent") %>%
  ggplot(aes(Mo17, sigma)) + geom_hex(bins=20) + facet_wrap(~type)+
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme

row2x <- ggplot_build(mv21)$layout$panel_ranges[[1]]$x.range
row2y <- ggplot_build(mv21)$layout$panel_ranges[[1]]$y.range

mv22 <- filter(mv.trend.df, type=="BNP") %>%
  ggplot(aes(Mo17, sigma)) + geom_hex(bins=20) + facet_wrap(~type)+
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme + xlim(row2x) + ylim(row2y)

mv31 <- filter(mv.trend.df, type=="independent") %>%
  ggplot(aes(flow_cell, sigma)) + geom_hex(bins=20) + facet_wrap(~type)+
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme

row3x <- ggplot_build(mv31)$layout$panel_ranges[[1]]$x.range
row3y <- ggplot_build(mv31)$layout$panel_ranges[[1]]$y.range

mv32 <- filter(mv.trend.df, type=="BNP") %>%
  ggplot(aes(flow_cell, sigma)) + geom_hex(bins=20) + facet_wrap(~type)+
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme + xlim(row3x) + ylim(row3y)


qqq2<- ggmatrix(list(mv11,mv12,
                     mv21,mv22,
                     mv31,mv32),3,2,xAxisLabels = c("independent","BNP"),
                yAxisLabels = c(paste(expression(paste(sigma," vs. intercept"))),
                                paste(expression(paste(sigma," vs. Mo17"))),
                                paste(expression(paste(sigma," vs. flow cell")))),
                labeller = label_parsed,
                showStrips = FALSE)#, legend = c(1,1)) + 
  theme(strip.text.y = element_text(size=12),
        strip.text.x = element_text(size=12))


cowplot::plot_grid(
  cowplot::plot_grid(mv1+my_theme+no.lgd,
                     mv2+my_theme+no.lgd,
                     mv3+my_theme+no.lgd,ncol=1),
  lgd, nrow=1, rel_widths = c(8,1))
  
ggsave("../figures_tables/mean-variance.pdf", width=8, height=11)

# New 9/13
ols12 <- select(ols, intercept, Mo17)
bnp12 <- select(bnp2, intercept, Mo17)
lmm12 <- select(lmm, intercept, Mo17)
source("../../chapter1/R/hist-diff.R")
p1 <- diffHexHist(bnp12,ols12, 30, .0008, "")+ggtitle("BNP-indep.")+no_panel+no_axes
p2 <- diffHexHist(lmm12,ols12, 30,.0008, "")+ggtitle("hierarch.-indep.")+theme(legend.position="none")+no_panel+no_axes
p3<- diffHexHist(bnp12,lmm12, 30,.0008, "")+ggtitle("BNP-hierarch.")+theme(legend.position = "none")+no_panel+no_axes
lgd <- cowplot::get_legend(p1)
p1 <- p1+theme(legend.position = "none") 
cowplot::plot_grid(cowplot::plot_grid(p1,p2,p3,ncol=1),
                   lgd, nrow=1)

ols13 <- select(ols, intercept, flow_cell)
bnp13 <- select(bnp2, intercept, flow_cell)
lmm13 <- select(lmm, intercept, flow_cell)
q1 <- diffHexHist(bnp13,ols13, 20, .0008, "")+ggtitle("BNP-indep.")+no_panel+no_axes
q2 <- diffHexHist(lmm13,ols13, 20,.0008, "")+ggtitle("hierarch.-indep.")+theme(legend.position="none")+no_panel+no_axes
q3<- diffHexHist(bnp13,lmm13, 20,.0008, "")+ggtitle("BNP-hierarch.")+theme(legend.position = "none")+no_panel+no_axes
lgd2 <- cowplot::get_legend(q1)
q1 <- q1+theme(legend.position = "none") 
cowplot::plot_grid(cowplot::plot_grid(q1,q2,q3,ncol=1),
                   lgd2, nrow=1, rel_widths = c(1,.17))

ols23 <- select(ols, Mo17, flow_cell)
bnp23 <- select(bnp2, Mo17, flow_cell)
lmm23 <- select(lmm, Mo17, flow_cell)
r1 <- diffHexHist(bnp23,ols23, 30, .0008, "")+ggtitle("BNP-indep.")+no_panel+no_axes
r2 <- diffHexHist(lmm23,ols23, 30,.0008, "")+ggtitle("hierarch.-indep.")+theme(legend.position="none")+no_panel+no_axes
r3<- diffHexHist(bnp23,lmm23, 30,.0008, "")+ggtitle("BNP-hierarch.")+theme(legend.position = "none")+no_panel+no_axes
lgd2 <- cowplot::get_legend(r1)
r1 <- r1+theme(legend.position = "none") 
cowplot::plot_grid(cowplot::plot_grid(r1,r2,r3,ncol=1),
                   lgd2, nrow=1, rel_widths = c(1,.17))


cowplot::plot_grid(p1,p2,p3,q1,q2,q3,r1,r2,r3, nrow=3)

ols14 <- data.frame(ols$intercept, sqrt(ind_est$sigma2-.1))
bnp14 <- select(bnp2, intercept, log_sigma) %>%
  mutate(sigma=exp(log_sigma)) %>%
  select(-log_sigma)
s1 <- diffHexHist(bnp14,ols14, 30, .0008, "")+ggtitle("BNP-OLS") +no_panel+no_axes+theme(legend.position="none")
s1

blank_panel <- ggally_blank()

p1 <- p1 + theme(legend.position="right")

ggmatrix(list(p1,p2,p3,q1,q2,q3,r1,r2,r3,s1, blank_panel, blank_panel),nrow = 4, ncol=3,yAxisLabels = c(paste(expression(paste(beta[2]," vs. ",beta[1]))),
                                                                           paste(expression(paste(beta[3],"  vs. ", beta[1]))),
                                                                           paste(expression(paste(beta[3],"  vs. ", beta[2]))),
                                                                           paste(expression(paste(sigma,"  vs. ", beta[1])))),
         xAxisLabels=c(paste(expression(paste("BNP", " vs. ", "independent"))),
                       paste(expression(paste("hierarchical", " vs. ", "independent"))),
                       paste(expression(paste("BNP", " vs. ", "hierarchical")))),
         labeller=label_parsed, legend=1)
ggsave("../figures_tables/difference_histograms.pdf", width=10, height=10)
