library(dplyr)
library(lme4)
library(ggplot2)

load("data/lmm_fit.RData")
load("data/ind_est.RData")
load("data/my_dat_wide.RData")

ols <- data.frame(t(ind_est$beta))
names(ols) <- c("intercept","halfdiff","flow_cell")
ols$GeneID <- my_dat_wide$GeneID
ols_long <- tidyr::gather(ols, key=genotype, value=est, -GeneID)
ols_long$type <- "ols"

lmm <- coef(lmm_fit)[[1]]
names(lmm)[c(1,3)] <- c("intercept","flow_cell")
lmm$GeneID <- my_dat_wide$GeneID
lmm_long <- tidyr::gather(lmm, key=genotype, value=est, -GeneID)
lmm_long$type <- "lmm"

bnp.fit <- readRDS("RDA-SD.rds")
bnp <- with(bnp.fit$summaries, data.frame(t(rbind(means_betas,log(means_sigmas)))))
names(bnp) <- c("intercept","halfdiff","flow_cell","log_sigma")
bnp$GeneID <- my_dat_wide$GeneID
bnp_long <- tidyr::gather(bnp, key=genotype, value=est, -GeneID)
bnp_long$type <- "BNP"

id <- sample(my_dat_wide$GeneID, size=40)

plot_df <- rbind(filter(ols_long),# GeneID %in% id),
                 filter(lmm_long),
                 filter(bnp_long,genotype!='log_sigma'))#, GeneID %in% id))



# ggplot(plot_df, aes(x=GeneID, y=est, color=type)) + geom_point(size=2) +
#   facet_wrap(~genotype, scales="free")
?getME

library(tidyr)
facet_names <- c( `lmm` = "LMM", `ols` = "OLS", `BNP` = "BNP",
                  `halfdiff` = "Half-Diff.",
                  `intercept` = "Intercept",
                  `flow_cell` = "Flow cell")

plot_df2 <- spread(plot_df, key=type, value=est)
p1 <- ggplot(plot_df2, aes(x=ols, y=lmm)) + geom_hex(bins=50) + 
  geom_abline(slope=1) + 
  facet_wrap(~genotype, scales="free", labeller=as_labeller(facet_names)) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))+
  xlab("OLS") + ylab("LMM")

p1.1 <- ggplot(plot_df2, aes(x=ols, y=BNP)) + geom_hex(bins=50) + 
  geom_abline(slope=1) + 
  facet_wrap(~genotype, scales="free", labeller=as_labeller(facet_names)) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))+
  xlab("OLS") + ylab("BNP")

p1.2 <- ggplot(plot_df2, aes(x=lmm, y=BNP)) + geom_hex(bins=50) + 
  geom_abline(slope=1) + 
  facet_wrap(~genotype, scales="free", labeller=as_labeller(facet_names)) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))+
  xlab("LMM") + ylab("BNP")

leg <- cowplot::get_legend(p1)
tm.noleg <- theme(legend.position = "none")
plt.grd1 <- cowplot::plot_grid(cowplot::plot_grid(p1+tm.noleg, p1.1+tm.noleg, p1.2+tm.noleg, ncol=1),leg,rel_widths = c(1,.15))

plot_df3 <- spread(plot_df, key=genotype, value=est)


p2 <- ggplot(plot_df3, aes(x=intercept, y=halfdiff)) + geom_hex() +
  facet_wrap(~type, labeller = as_labeller(facet_names))+
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))



plot_df_ws <- cbind(ols, log_sigma = .5*log(ind_est$sigma2)) %>%
  mutate(GeneID = my_dat_wide$GeneID)

library(GGally)

### Need to set common scales manually
my_theme <- theme_bw(base_size=12) + theme(panel.grid = element_blank())
q11 <- filter(plot_df3, type=="ols") %>%
  ggplot(aes(x=intercept, y=halfdiff)) + geom_hex() +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col1xlim <- ggplot_build(q11)$layout$panel_ranges[[1]]$x.range
col1ylim <- ggplot_build(q11)$layout$panel_ranges[[1]]$y.range


q12 <- filter(plot_df3, type=="ols") %>%
  ggplot(aes(x=intercept, y=flow_cell)) + geom_hex() +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col2xlim <- ggplot_build(q12)$layout$panel_ranges[[1]]$x.range
col2ylim <- ggplot_build(q12)$layout$panel_ranges[[1]]$y.range


q13 <- filter(plot_df3, type=="ols") %>%
  ggplot(aes(x=halfdiff, y=flow_cell)) + geom_hex() +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col3xlim <- ggplot_build(q13)$layout$panel_ranges[[1]]$x.range
col3ylim <- ggplot_build(q13)$layout$panel_ranges[[1]]$y.range
# q13 <- q13 + ylim(col3ylim) + xlim(col3xlim)

q21 <- filter(plot_df3, type=="lmm") %>%
  ggplot(aes(x=intercept, y=halfdiff)) + geom_hex() +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme + xlim(col1xlim) + ylim(col1ylim)

q22 <- filter(plot_df3, type=="lmm") %>%
  ggplot(aes(x=intercept, y=flow_cell)) + geom_hex() + 
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme + xlim(col2xlim) + ylim(col2ylim)
  # #geom_violin(color="darkblue", alpha=0)+
  # geom_vline(xintercept=log(summary(lmm_fit)$sigma), linetype=2, color="darkblue")+
  # my_theme + coord_flip() + ylim(col2xlim) + xlim(col2ylim)

q23 <- filter(plot_df3, type=="lmm") %>%
  ggplot(aes(x=halfdiff, y=flow_cell)) + geom_hex() + 
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme + xlim(col3xlim) + ylim(col3ylim)
# %>%
#   ggplot(aes(x=log(summary(lmm_fit)$sigma), y=halfdiff)) + geom_violin(color="darkblue", alpha=0)+
#   geom_vline(xintercept=log(summary(lmm_fit)$sigma), linetype=2, color="darkblue")+
#   my_theme + coord_flip() + ylim(col3xlim) + xlim(col3ylim)

# qnull <- ggplot() + theme_void()
q31 <- bnp %>%
  ggplot(aes(x=intercept, y=halfdiff)) + geom_hex() +
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme + xlim(col1xlim) + ylim(col1ylim)

q32 <- bnp %>%
  ggplot(aes(x=intercept, y=flow_cell)) + geom_hex() +
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme + xlim(col2xlim) + ylim(col2ylim)

q33 <- bnp %>%
  ggplot(aes(x=halfdiff, y=flow_cell)) + geom_hex() +
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  my_theme + xlim(col3xlim) + ylim(col3ylim)

qqq <- ggmatrix(list(q11,q12,q13,q21,q22,q23,q31,q32,q33),3,3,
                xAxisLabels = c("paste('Half-Diff. vs Intercept')",
                                "paste('Flow Cell vs Intercept')",
                                "paste('Flow Cell vs Half-Diff.')"),
                yAxisLabels = c("OLS", "LMM", "BNP"),
                labeller = label_parsed,
                legend=1
) +theme(strip.text.y = element_text(size=12),
         strip.text.x = element_text(size=12))

#volcanoplot
volc.df <- data.frame(effect = bnp.fit$summaries$means_betas[3,],
                      probs = bnp.fit$summaries$probs[1,])
ggplot(volc.df, aes(effect,probs)) + geom_hex(bins=50) +
  scale_fill_continuous(trans="log",breaks=c(1,10,100,1000),low="white", high = "darkblue")+
  geom_vline(xintercept=0, linetype=2)
