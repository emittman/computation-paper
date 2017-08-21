library(Paschold2012)
library(dplyr)
library(tidyr)
source("getNormFactors.R")

#differential expression, B73 v. B73xMo17
my_dat_wide <- Paschold2012 %>%
  mutate(genotype_replicate = paste(genotype,replicate,sep="_")) %>%
  select(GeneID, genotype_replicate, total) %>%
  tidyr::spread(genotype_replicate, total) %>%
  select(1:9)

#remove all zero genes
zero_id <- which(apply(data.matrix(my_dat_wide[,-1]), 1, function(x) all(x==0)))
my_dat_wide <- my_dat_wide[-zero_id,]

#transform counts to log (+1) scale, normalize
my_dat_wide[,2:9] <- normalizeData(data.matrix(my_dat_wide[,2:9]),
                                   group=rep(1:2, each=4), trans.to.log = TRUE)

save(my_dat_wide, file="data/my_dat_wide.RData")

#convert to long format, add indicators for genotype and flowcell
my_dat <- gather(my_dat_wide, key=genotype_replicate, value=total, -GeneID) %>%
  extract(col=genotype_replicate, into=c("genotype","replicate"),
                     regex = "([[:alnum:]]+)_([[:alnum:]]+)", convert=TRUE) %>%
  mutate(halfdiff = ifelse(genotype=="B73",-1,1),
         flow_cell = factor(ifelse(replicate %in% c(1,2), 1, 2))) %>%
  arrange(GeneID)

save(my_dat, file="data/my_dat.RData")

source("../../../cuda_rpackage/R/data.R")

X <- filter(my_dat, GeneID == my_dat$GeneID[1]) %>% head(8) %>%
  model.matrix(~halfdiff+flow_cell, data=.)

y <- data.matrix(my_dat_wide[,2:9])

cuda_dat <- formatData(counts = y, groups = rep(1:2, each=4), X = X, voom = FALSE, transform_y=identity)
ind_est <- indEstimates(cuda_dat)

save(ind_est, file="data/ind_est.RData")
save(cuda_dat, file="data/cuda_dat.Rdata")
