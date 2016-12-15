#######################################
#
#	Codes from EMMA

#
#	log restricted log-likelihood function
	
emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
	nq <- length(etas)
	delta <-  exp(logdelta)
	return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}


emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

## Core functions of EPACTS to perform association

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

single.q.emmaxR <- function() {
  ## fit the null model adjusting to covariates
  cname <- c("STAT","BETA","SEBETA","R2")
  print("foo")
  eig <- .Call("readEigenWithIDs",eigf);
  print("bar")

  ## read REML output
  remls <- as.matrix(read.table(remlf)[,2])
  delta <- remls[1]
  #va <- remls[2]
  #ve <- remls[3]

  ## T <- evecs * sqrt(1/diag(evals+delta))
  n <- nrow(eig$evecs)
  q <- ncol(eig$evecs)
  emxT <- matrix( 1/sqrt(eig$evals+delta), q, n ) * t(eig$evecs)  ## not matrix multiplication

  ## run score test
  phenoT <- emxT %*% pheno  ## q * 1 matrix
  genoT <- emxT %*% t(genos)    ## q * m matrix

  ## run simple linear regression
  sy <- sum(phenoT)
  syy <- sum(phenoT*phenoT)
  sx <- colSums(genoT)
  sxx <- colSums(genoT*genoT)
  sxy <- crossprod(phenoT, genoT)

  beta <- ( (q+1)*sxy - sx * sy)/ ( (q+1)*sxx - sx*sx )
  varE <- 1/(q+1)/(q-1)*( (q+1)*syy - sy*sy - beta*beta*((q+1)*sxx-sx*sx))
  sebeta <- sqrt((q+1)*varE/((q+1)*sxx-sx*sx))
  r <- ((q+1)*sxy-sx*sy)/sqrt(((q+1)*sxx-sx*sx)*((q+1)*syy-sy*sy))
  t <- r * sqrt((q-1)/(1- r*r + 1e-10))
  #print(n)
  #print(q)
  #stop()
  return(list(p=pt(abs(t),q-1,lower.tail=F)*2,add=t(rbind(t,beta,sebeta,r*r)),cname=cname))
}

