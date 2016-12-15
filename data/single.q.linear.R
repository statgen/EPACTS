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

## single.q.linear() : Single Variant Linear Regression
## KEY FEATURES :
##  1. Should be faster than lm() implementation
##  2. Not exact : regress out the covariates first
## TRAITS  : QUANTITATIVE (GAUSSIAN)
## RETURNS : PVALUE, TSTAT, BETA, SEBETA, R2
## MISSING VALUE : IGNORED
single.q.linear <- function() {
  cname <- c("BETA","SEBETA","TSTAT","R2")
  res <- lm(pheno ~ cov - 1)$residual
  vAC <- AC[vids]
  vNS <- NS[vids]
  df0 <- vNS-1
  #print(pheno)
  #print(cov)
  #print(genos[1,])
  #stop()  
  
  genos <- genos - matrix(vAC/vNS,nrow(genos),n) ## make mean zero
  s2x <- sqrt(rowSums(genos * genos, na.rm=T)/df0)
  if ( sum(n-vNS) > 0 ) { ## if missing genotype exists
    #print(length(res))
    #print(dim(genos))
    phenos <- matrix(res,nrow(genos),n,byrow=T)
    phenos[is.na(genos)] <- NA
    mean.phenos <- rowMeans(phenos,na.rm=T)
    s2y <- sqrt(rowSums(phenos * phenos, na.rm=T)/df0)
    sxy <- rowSums(genos * phenos, na.rm=T)/df0
  }
  else {
    s2y <- sqrt(sum(res * res)/(n-1))
    sxy <- tcrossprod(res, genos)/df0
  }
  rxy <- sxy / ( s2y * s2x )
  b <- rxy * ( s2y / s2x )
  se.b <- s2y / s2x * sqrt ( (1-rxy*rxy+1e-20) / (df0-1) )
  t <- b/se.b
  return(list(p=pt(abs(t),vNS-2,lower.tail=F)*2,
              add=matrix(c(b,se.b,t,rxy*rxy),nrow(genos),4),
              cname=cname))
}

