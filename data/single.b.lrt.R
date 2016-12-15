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

## single.b.lrt:
## Likelihood ratio test
## Use built-in glm logistic regression to perform association
## By:  Clement Ma
##
## KEY FEATURES :SIMPLE, BUT MAY BE SLOW
## TRAITS  : BINARY
## RETURNS : PVALUE, BETA, SEBETA, ZSTAT
## MISSING VALUES : IGNORED
single.b.lrt <- function() {
  #cname <- NULL
  m <- nrow(genos)
  p <- rep(NA,m)
  cname <- c("NULL_DEVIANCE","FULL_DEVIANCE","N.CASE","N.CTRL","AF.CASE","AF.CTRL");
  add <- matrix(NA,m,6)
  #add <- matrix(NA,m,3) ## BETA, SEBETA, TSTAT
  
  # Model without genotype
  r.nogeno <- glm(pheno~cov-1, family="binomial")
  if ( m > 0 ) {
    for(i in 1:m) {
      # add back genotype
      r <- update(r.nogeno, ~ . + t(genos[i,,drop=FALSE]))
      if ( ( r$converged ) && ( ! r$boundary ) ) {
        #p[i] <- anova(r, r.nogeno, test="Chisq")[2,"Pr(>Chi)"]
        p[i] <- pchisq(r.nogeno$deviance-r$deviance, df=1, lower.tail=F)
        add[i,1:2] <- c(r.nogeno$deviance,r$deviance)
      }
      else {
        p[i] <- NA
      }
    }
  }
  ina <- is.na(genos)
  n.cases <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 1) * (1-ina))
  n.ctrls <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 0) * (1-ina))
  ac.cases <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 1) * genos)
  ac.ctrls <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 0) * genos)    
  
  add[,3] <- n.cases
  add[,4] <- n.ctrls
  add[,5] <- ac.cases/(2*n.cases)
  add[,6] <- ac.ctrls/(2*n.ctrls)
  return(list(p=p,add=add,cname=cname))
}
