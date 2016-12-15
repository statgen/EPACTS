## Core functions of pugwash to perform association

##################################################################
## SINGLE VARIANT TEST
## INPUT VARIABLES:
##   n        : total # of individuals
##   NS       : number of called samples
##   AC       : allele count
##   MAF      : minor allele frequency
##   vids     : indices from 1:n after AF/AC threshold
##   genos    : genotype matrix (after AF/AC threshold)
## EXPECTED OUTPUT : list(p, addcols, addnames) for each genos row
##   p        : p-value
##   addcols  : additional columns (with proper column names)
##################################################################

## multi.q.linear() : Single Variant Linear Regression
## KEY FEATURES :
##  1. Should be faster than lm() implementation
##  2. Not exact : regress out the covariates first
## TRAITS  : QUANTITATIVE (GAUSSIAN)
## RETURNS : PVALUE, TSTAT, BETA, SEBETA, R2
## MISSING VALUE : IGNORED
multi.q.linear <- function() {  ## note that phenotypes are already regressed out
  ## make mean zero
  genos.mean <- matrix(rowMeans(genos,na.rm=T),nrow(genos),ncol(genos),byrow=F)
  genos[is.na(genos)] <- genos.mean[is.na(genos)] ## Mean imputation

  vNS <- NS[vids]
  df0 <- vNS-1

  ## phenos are n * g matrix, genos are m * n matrix
  T <- cor(t(genos),phenos)
  s2y <- sqrt(colSums(phenos * phenos)/(n-1))
  s2x <- sqrt(rowSums(genos * genos, na.rm=T)/df0)
  B <- matrix(s2y,nrow(T),ncol(T),byrow=T) / matrix(s2x,nrow(T),ncol(T),byrow=F) * T
  T <- T / sqrt((1-T*T+1e-20)/(df0-1))

  return(list(p=matrix(pt(abs(T),df0-1,lower.tail=F)*2,nrow(B),ncol(B)), b=B))
}

