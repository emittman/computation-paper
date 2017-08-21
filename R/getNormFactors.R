library(edgeR)

getNormFactors <- function(counts, group){
  edgeR_dat <- edgeR::DGEList(counts=counts, group=group)
  edgeR_dat <- calcNormFactors(edgeR_dat)
  edgeR_dat$samples$norm.factors
}

normalizeData <- function(counts, group, trans.to.log=T){
  normFactors <- getNormFactors(counts, group)
  print(normFactors)
  if(trans.to.log){
    return(sapply(1:ncol(counts), function(i) log(counts[,i]+1) + log(normFactors[i])))
  } else{
    sapply(1:ncol(counts), function(i) counts[,i] * normFactors[i])
  }
}
