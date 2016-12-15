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

## single.q.reverse() : Reverse logistic regression
## KEY FEATURES : ASSUMPTION OF DOMINANT MODEL FOR MINOR ALLELE
##                COVARIATES ARE ADJUSTED TOGETHER
## TRAITS  : QUANTITATIVE (GAUSSIAN)
## RETURNS : PVALUE, BETA, SEBETA, ZSTAT
## MISSING VALUE : IGNORED
single.q.reverse <- function() {
  cname <- c("BETA","SEBETA","ZSTAT")
  m <- nrow(genos)
  vflip <- which(AC[vids]>NS[vids])
  genos[vflip,,drop=FALSE] = 2-genos[vflip,,drop=FALSE] # flip  
  p <- rep(NA,m)
  add <- matrix(NA,m,3)
  for(i in 1:m) {
    g <- (t(genos[i,,drop=FALSE]) > 0)
    r <- glm(g~pheno+cov-1,family=binomial)
    if ( ( r$converged ) && ( ! r$boundary ) ) {
      p[i] <- summary(r)$coefficients[1,4]
      add[i,] <- summary(r)$coefficients[1,1:3]
    }
  }
  return(list(p=p,
              add=add,
              cname=cname))
}
