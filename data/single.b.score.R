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

## single.b.logit.score() :
##   Single Variant Score Test for Binary Traits
## # Modified from Clement Ma's software
## TRAITS  : BINARY ONLY
## RETURNS : PVALUE, SCORES
## MISSING VALUE : IMPUTATION BY MEAN
single.b.score <- function() {
  ## fit the null model adjusting to covariates
  cname <- c("SCORE")
  #cname <- c("SCORE","N.CASE","N.CTRL","AF.CASE","AF.CTRL")  
  reml <- glm(pheno ~ cov - 1, family=binomial, x=T)  
  ## calculate the weights for each observation
  v <- exp(reml$linear.predictors) / (1 + exp(reml$linear.predictors))^2
  vx <- v * reml$x
  v2 <- t(vx) %*% reml$x
  iv2 <- solve(v2)
  iv2c <- chol(iv2) ## make t(iv2c) %*% iv2c = iv2
  x.res <- pheno - reml$fitted  # n * 1 matrix

  ## resolve missing genotype by mean imputation
  ina <- is.na(genos)
  if ( length(vids) > 0 ) {
    genos[ina] <- matrix(AC[vids]/NS[vids],nrow(genos),ncol(genos))[ina]
  }
  #n.cases <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 1) * (1-ina))
  #n.ctrls <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 0) * (1-ina))
  #ac.cases <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 1) * genos * (1-ina))
  #ac.ctrls <- rowSums((matrix(pheno,nrow(genos),ncol(genos),byrow=T) == 0) * genos * (1-ina))

  ## Multiply weights for the genotypes
  U <- tcrossprod(t(x.res),genos)  # U is m * 1 matrix
  Vl <- tcrossprod(v,genos*genos)
  Vr <- genos %*% tcrossprod(vx,iv2c)
  V.s <- sqrt(Vl - rowSums(Vr*Vr) + 1e-10) ## to avoid negative variance estimate
  #print(min(Vl-rowSums(Vr*Vr)))

  ## compute pvalues and scores and return them
  T <- U/V.s
  return(list(p=pnorm(abs(T),lower.tail=F)*2,add=cbind(matrix(T,length(T),1)),cname=cname))
  #return(list(p=pnorm(abs(T),lower.tail=F)*2,add=cbind(matrix(T,length(T),1),n.cases,n.ctrls,ac.cases/(2*n.cases),ac.ctrls/(2*n.ctrls)),cname=cname))
}
