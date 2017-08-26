#Good agreement between two methods

setwd("R")
sd <- readRDS("RDA-SD.rds")
sb <- readRDS("RDA-SB.rds")
library(ggplot2)
library(dplyr)

beta_sd <- data.frame(with(sd, t(summaries$means_betas))) %>%
  mutate(gene = 1:nrow(beta_sd), type="symmetric-dirichlet")
names(beta_sd)[1:3] <- c("intercept", "half-diff", "flow-cell")
beta_sb <- data.frame(with(sb, t(summaries$means_betas))) %>%
  mutate(gene = 1:nrow(beta_sd), type="stick-breaking")
names(beta_sb)[1:3] <- c("intercept", "half-diff", "flow-cell")

melt_df <- rbind(beta_sd,beta_sb) %>%
  tidyr::gather(key=effect, value=value, 1:3)


melt_df %>% filter(effect=="intercept") %>%
  tidyr::spread(key=type, value=value) %>% 
  ggplot(aes_string(x="`stick-breaking`",y="`symmetric-dirichlet`")) + 
  geom_hex(bins=80) + geom_abline(slope=1)

melt_df %>% filter(effect=="half-diff") %>%
  tidyr::spread(key=type, value=value) %>% 
  ggplot(aes_string(x="`stick-breaking`",y="`symmetric-dirichlet`")) + 
  geom_hex(bins=80) + geom_abline(slope=1)

melt_df %>% filter(effect=="flow-cell") %>%
  tidyr::spread(key=type, value=value) %>% 
  ggplot(aes_string(x="`stick-breaking`",y="`symmetric-dirichlet`")) + 
  geom_hex(bins=80) + geom_abline(slope=1)
