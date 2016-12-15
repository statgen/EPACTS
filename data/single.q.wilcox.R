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

## single.wilcox() : Wilcoxon Rank-sum Test for Quantitative Traits
## KEY FEATURES : ASSUMPTION OF DOMINANT MODEL FOR MINOR ALLELE
##                NON-PARAMETRIC, MAY BE SLOW
##                COVARIATES ARE ADJUSTED BEFORE RUNNING
## TRAITS  : QUANTITATIVE 
## RETURNS : PVALUE, STAT
## MISSING VALUE : IGNORED
single.q.wilcox <- function() {
  cname <- c("STAT")  
  res <- lm(pheno ~ cov - 1)$residual
  m <- nrow(genos)
  vflip <- which(AC[vids]>NS[vids])
  genos[vflip,] = 2-genos[vflip,] # flip
  p <- rep(NA,m)
  T <- rep(NA,m)
  for(i in 1:m) {
    g <- (genos[i,] > 0) # dominant model
    r <- wilcox.test(res~g)
    p[i] <- r$p.value
    T[i] <- r$statistic
  }
  return(list(p=p,
              add=matrix(T,m,1),
              cname=cname))
}
