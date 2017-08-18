library(dplyr)
library(lme4)
library(ggplot2)

lmm_fit <- readRDS("data/lmm_fit.rds")
est <- readRDS("data/est.rds")

load("diff_expr.Rdata")

ols <- data.frame(t(est$beta))
names(ols) <- c("intercept","halfdiff")
ols$GeneID <- my_dat_wide$GeneID
ols_long <- tidyr::gather(ols, key=genotype, value=est, -GeneID)
ols_long$type <- "ols"

lmm <- coef(lmm_fit)[[1]]
names(lmm)[1] <- "intercept"
lmm$GeneID <- my_dat_wide$GeneID
lmm_long <- tidyr::gather(lmm, key=genotype, value=est, -GeneID)
lmm_long$type <- "lmm"

id <- sample(wide$GeneID, size=40)

plot_df <- rbind(filter(ols_long),# GeneID %in% id),
                 filter(lmm_long))#, GeneID %in% id))



library(ggplot2)
ggplot(plot_df, aes(x=GeneID, y=est, color=type)) + geom_point(size=2) +
  facet_wrap(~genotype, scales="free")
?getME

library(tidyr)
facet_names <- c( `lmm` = "LMM", `ols` = "OLS",
                  `halfdiff` = "Half-Diff.",
                  `intercept` = "Intercept")

plot_df2 <- spread(plot_df, key=type, value=est)
p1 <- ggplot(plot_df2, aes(x=ols, y=lmm)) + geom_hex() + 
  geom_abline(slope=1) + 
  facet_wrap(~genotype, scales="free", labeller=as_labeller(facet_names)) +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))+
  xlab("OLS") + ylab("LMM")

plot_df3 <- spread(plot_df, key=genotype, value=est)


p2 <- ggplot(plot_df3, aes(x=intercept, y=halfdiff)) + geom_hex() +
  facet_wrap(~type, labeller = as_labeller(facet_names))+
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  theme_bw(base_size=12)+theme(strip.text.y = element_text(size=12),
                               strip.text.x = element_text(size=12))



plot_df_ws <- cbind(ols, log_sigma = .5*log(est$sigma2)) %>%
  mutate(GeneID = my_dat_wide$GeneID)

library(GGally)

### Need to set common scales manually
my_theme <- theme_bw(base_size=12)
q11 <- filter(plot_df3, type=="ols") %>%
  ggplot(aes(x=intercept, y=halfdiff)) + geom_hex() +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col1xlim <- ggplot_build(q11)$layout$panel_ranges[[1]]$x.range
col1ylim <- ggplot_build(q11)$layout$panel_ranges[[1]]$y.range


q12 <- ggplot(plot_df_ws, aes(x=intercept, y=log_sigma)) + geom_hex() +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col2xlim <- ggplot_build(q12)$layout$panel_ranges[[1]]$x.range
col2ylim <- ggplot_build(q12)$layout$panel_ranges[[1]]$y.range

q13 <- ggplot(plot_df_ws, aes(x=halfdiff, y=log_sigma)) + geom_hex() +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme

col3xlim <- ggplot_build(q13)$layout$panel_ranges[[1]]$x.range
col3ylim <- ggplot_build(q13)$layout$panel_ranges[[1]]$y.range

q21 <- filter(plot_df3, type=="lmm") %>%
  ggplot(aes(x=intercept, y=halfdiff)) + geom_hex() +
  scale_fill_continuous(trans="log", breaks=c(1, 10, 100, 1000), low="white", high="darkblue")+
  my_theme + xlim(col1xlim) + ylim(col1ylim)

q22 <- filter(plot_df3, type=="lmm") %>%
  ggplot(aes(x=log(summary(lmm_fit)$sigma), y=intercept)) + geom_violin(color="darkblue", alpha=0)+
  geom_vline(xintercept=log(summary(lmm_fit)$sigma), linetype=2, color="darkblue")+
  my_theme + coord_flip() + ylim(col2xlim) + xlim(col2ylim)

q23 <- filter(plot_df3, type=="lmm") %>%
  ggplot(aes(x=log(summary(lmm_fit)$sigma), y=halfdiff)) + geom_violin(color="darkblue", alpha=0)+
  geom_vline(xintercept=log(summary(lmm_fit)$sigma), linetype=2, color="darkblue")+
  my_theme + coord_flip() + ylim(col3xlim) + xlim(col3ylim)

# qnull <- ggplot() + theme_void()

library(GGally)

qqq <- ggmatrix(list(q11,q12,q13,q21,q22,q23),2,3,
                xAxisLabels = c("paste('Half-Diff. vs Intercept')",
                                "paste(log(sigma),' vs Intercept')",
                                "paste(log(sigma),' vs Half-Diff.')"),
                yAxisLabels = c("OLS", "LMM"),
                labeller = label_parsed,
                legend=1
) +theme(strip.text.y = element_text(size=12),
         strip.text.x = element_text(size=12))


