library(Paschold2012)
library(dplyr)
my_dat <- Paschold2012::Paschold2012

#differential expression, B73 v. B73xMo17
my_dat_wide <- Paschold2012 %>%
  mutate(genotype_replicate = paste(genotype,replicate,sep="_")) %>%
  select(GeneID, genotype_replicate, total) %>%
  tidyr::spread(genotype_replicate, total) %>%
  select(1:9)

#remove all zero genes
zero_id <- which(apply(data.matrix(my_dat_wide[,-1]), 1, function(x) all(x==0)))

my_dat <- my_dat[which(!(my_dat$GeneID %in% my_dat_wide$GeneID[zero_id])),]
my_dat_wide <- my_dat_wide[-zero_id,]

my_dat <- mutate(my_dat,
                 halfdiff = ifelse(genotype=="B73",-1,1))

save(my_dat, my_dat_wide, file="diff_expr.Rdata")