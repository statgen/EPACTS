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

## single.lm() : Use built-in lm() function to perform association
## KEY FEATURES : SIMPLE, BUT MAY BE SLOW
##                GOOD SNIPPLET TO START A NEW FUNCTION
## TRAITS  : QUANTITATIVE
## RETURNS : PVALUE, BETA, SEBETA, TSTAT
## MISSING VALUES : IGNORED
single.q.lm <- function() {
  cname <- c("BETA","SEBETA","TSTAT")
  m <- nrow(genos)
  p <- rep(NA,m)
  add <- matrix(NA,m,3) ## BETA, SEBETA, TSTAT
  if ( m > 0 ) {
    for(i in 1:m) {
      print(dim(genos))
      r <- summary(lm(pheno~t(genos[i,,drop=FALSE])+cov-1))$coefficients[1,]
      p[i] <- r[4]
      add[i,] <- r[1:3]
    }
  }
  return(list(p=p,add=add,cname=cname))
}
