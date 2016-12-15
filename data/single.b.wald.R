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

## single.b.logit.wald :
## Use built-in wald logistic regression to perform association
## KEY FEATURES :SIMPLE, BUT MAY BE SLOW
## TRAITS  : BINARY
## RETURNS : PVALUE, BETA, SEBETA, ZSTAT
## MISSING VALUES : IGNORED
single.b.wald <- function() {
  #cname <- c("BETA","SEBETA","ZSTAT","N.CASE","N.CTRL","AF.CASE","AF.CTRL")
  cname <- c("BETA","SEBETA","ZSTAT")
  m <- nrow(genos)
  p <- rep(NA,m)
  add <- matrix(NA,m,length(cname)) ## BETA, SEBETA, TSTAT
  if ( m > 0 ) {
    for(i in 1:m) {
      r <- glm(pheno~t(genos[i,,drop=FALSE])+cov-1,family="binomial")
      if ( ( r$converged ) && ( ! r$boundary ) ) {
        s <- summary(r)$coefficients[1,]
        p[i] <- s[4]
        add[i,1:3] <- s[1:3]
      }
      else {
        p[i] <- NA
        add[i,1:3] <- rep(NA,3)
      }
    }
  }
  ## resolve missing genotype by mean imputation
  ##ina <- is.na(genos)
  ##if ( length(vids) > 0 ) {
    ##genos[ina] <- matrix(AC[vids]/NS[vids],nrow(genos),ncol(genos))[ina]
  ##}
  ##n.cases <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 1) * (1-ina))
  ##n.ctrls <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 0) * (1-ina))
  ##ac.cases <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 1) * genos * (1-ina))
  ##ac.ctrls <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 0) * genos * (1-ina))
  
  ##add[,4] <- n.cases
  ##add[,5] <- n.ctrls
  ##add[,6] <- ac.cases/(2*n.cases)
  ##add[,7] <- ac.ctrls/(2*n.ctrls)
  return(list(p=p,add=add,cname=cname))
}
