## check the installation of SKAT package
if ( !require(mmSKAT,lib.loc=paste(bindir,"/lib/",sep="") ) ) {
  stop("Cannot find mmSKAT package");
}

library(mmSKAT)


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

#######################
#
# Changed by SLEE

SKAT.emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  XX1<-X %*% (solve(crossprod(X,X)))
  K1<-t(X) %*% K
  K2<-K %*% X
  
  #Mat<-K - XX1 %*% K1 - K2 %*% t(XX1) + XX1 %*% ((K1 %*% X) %*% t(XX1)) - XX1 %*% t(X)
  Mat<-K - K2 %*% t(XX1) - XX1 %*% (K1  +  ((K1 %*% X) %*% t(XX1)) -  t(X))
  
  diag(Mat) = diag(Mat) +1
  eig <- eigen(Mat,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

#
#	X : covariates
#	Z = NULL
#	
SKAT_NULL_emmaX <- function(formula, data=NULL, K=NULL, Kin.File=NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
  
        # check missing 
        obj1<-model.frame(formula,na.action = na.pass,data)
	obj2<-model.frame(formula,na.action = na.pass,data)

	n1<-dim(obj2)[1]
	n<-dim(obj1)[1]
	id_include<-SKAT:::SKAT_Null_Model_Get_Includes(obj1,obj2)
	X<-model.matrix(formula,data=data)
	y<-model.response(obj1)
	q <- ncol(X)
	# n and n1 are opposite (compare to SKAT)
	if(n1 - n > 0){
		MSG<-sprintf("%d  samples have either missing phenotype or missing covariates. They are excluded from the analysis!",n1 - n)
		warning(MSG,call.=FALSE)
	}

	stopifnot(nrow(X) == n)
	
	########################################################
	# Read kinship 
	
	if(is.null(K) && is.null(Kin.File)){
		stop("both K and Kin.File are NULL!")
	} else if(is.null(K)){
		if(!file.exists(Kin.File)){
			msg<-sprintf("File %s is not exist!", Kin.File)
			stop(msg)
		}
		cat("Read ", Kin.File, "\n")
		K = as.matrix(read.delim(Kin.File, header=FALSE, colClasses="numeric", sep = "\t"))
		t <- nrow(K)
		cat("Read complete. ", Kin.File, " has ", t, " rows! \n")
	}
	
	if(class(K) != "matrix"){
		stop("K is not a matrix!")
	}
	
	if(n1 - n > 0){
		K<-K[id_include, id_include]
	}
	t <- nrow(K)
	stopifnot(ncol(K) == t)
		
	#######################################################
	# Estimate parameters
		
    eig.R <- emma.eigen.R.wo.Z(K,X)
	
    etas <- crossprod(eig.R$vectors,y)
  	logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      	optlogdelta <- append(optlogdelta, llim)
      	optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      	optlogdelta <- append(optlogdelta, ulim)
      	optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }

    for( i in 1:(m-1) ) {
        if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ){
          	r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
          	optlogdelta <- append(optlogdelta, r$root)
          	optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
              }
      }

	#############################################################################
  	# variance term : maxva K + maxve I
  	# eig.R$vectors $*$ diag(maxva * eig.R$values  + maxdelta) %*% t(eig.R$vectors)
  	
  	
	maxdelta <- exp(optlogdelta[which.max(optLL)])
  	maxLL <- max(optLL)
	maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)    # additive effect, genetic

  	maxve <- maxva*maxdelta	# noise
  	 
	#############################################################################
  	# Get V^-1 (y - X \beta)
  	# = V-1 y - V-1 X (X'V-1X)-1 X'V-1 y
  	# = P y
  	# P = V-1 -  V-1 X (X'V-1X)-1 X'V-1
  
    #lambda_inv<-1/(maxva*eig.R$values+maxve)
  	#X.etas <- t(eig.R$vectors) %*% X 
    #XVX<- t(X.etas) %*% (X.etas * lambda_inv)
  	#XVX_inv<-solve(XVX)		
  	#V_inv<- eig.R$vectors %*% (t(eig.R$vectors) * lambda_inv )
  	#XV_inv = t(X) %*% V_inv
  	#P = V_inv -  t(XV_inv) %*% (XVX_inv %*% XV_inv)
  	#res = P %*% y
 	
 	va=maxva
 	ve=maxve 	
  	
  	V = va * K
  	diag(V) = diag(V) + ve 
    V_inv <- solve(V)
  
  	XVX = t(X) %*% (V_inv %*% X)
  	XVX_inv<-solve(XVX)
  	
  	
  	XV_inv = t(X) %*% V_inv
  	P = V_inv -  t(XV_inv) %*% (XVX_inv %*% XV_inv)
  	res = P %*% y
  	
  	
  	re<-list( LL=maxLL, va=va, ve=ve, P=P, res=res, id_include=id_include)
  	class(re)<-"SKAT_NULL_Model_EMMAX"
  	return (re)
}

SKAT_emmaX = function( Z, obj, kernel= "linear.weighted", weights.beta=c(1,25), weights=NULL, method="davies", impute.method="fixed", r.corr=0, is_check_genotype=TRUE, is_dosage = FALSE, missing_cutoff=0.15, SetID = NULL){
  if(class(obj) != "SKAT_NULL_Model_EMMAX"){
    stop("ERROR: obj is not an returned object from SKAT.emmaX.null")
  }
  
  m = ncol(Z)
  n = nrow(Z)
#####################################
# Check genotypes and parameters
	out.z<-SKAT:::SKAT_MAIN_Check_Z(Z, n, obj$id_include, SetID, weights, weights.beta, impute.method, is_check_genotype, is_dosage, missing_cutoff)
	if(out.z$return ==1){
		out.z$param$n.marker<-m
		return(out.z)
	}
	Z = out.z$Z.test
	weights = out.z$weights
  	
  	# Weighted Linear Kernel 
  	if (kernel == "linear.weighted") {
    	Z = t(t(Z) * (weights))
  	}
  
  	if(r.corr == 1){
  		Z<-cbind(rowSums(Z))
  	} else if(r.corr > 0){

   		p.m<-dim(Z)[2]	
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		L<-chol(R.M,pivot=TRUE)
		Z<- Z %*% t(L) 
  	}

  	# get Q
  	Q.Temp = t(obj$res)%*%Z
  	Q = Q.Temp %*% t(Q.Temp)/2
  	Q.res=NULL

	# Get Z' P0 Z
  	W.1= t(Z) %*% (obj$P %*% Z) # t(Z) P0 Z

 	if( method == "liu.mod" ){
		out<-SKAT:::Get_Liu_PVal.MOD(Q, W.1, Q.res)    

  	} else if( method == "davies" ){

		out<-SKAT:::Get_Davies_PVal(Q, W.1, Q.res)    

  	} else {
		stop("Invalid Method!")
  	}

  
  re<-list(p.value = out$p.value, Test.Type = method, Q = Q, param=out$param ) 

  re$param$n.marker<-m
  re$param$n.marker.test<-ncol(out.z$Z.test)
	 
  return(re)
}

#
# x is either y or SKAT_NULL_Model 
#
SKAT_emmaX.SSD.OneSet = function(SSD.INFO, SetID, obj, ...){
	
	id1<-which(SSD.INFO$SetInfo$SetID == SetID)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set id [%s] from SSD!", SetID)
		stop(MSG)
	}	
	Set_Index<-SSD.INFO$SetInfo$SetIndex[id1]

	Z<-SKAT:::Get_Genotypes_SSD(SSD.INFO, Set_Index)
	re<-SKAT_emmaX(Z, obj, ...)
	
	return(re)
}

#
# x is either y or SKAT_NULL_Model 
#
SKAT_emmaX.SSD.OneSet_SetIndex = function(SSD.INFO, SetIndex, obj, ...){
	
	id1<-which(SSD.INFO$SetInfo$SetIndex == SetIndex)
	if(length(id1) == 0){
		MSG<-sprintf("Error: cannot find set index [%d] from SSD!", SetIndex)
		stop(MSG)
	}	
	SetID<-SSD.INFO$SetInfo$SetID[id1]


	Z<-SKAT:::Get_Genotypes_SSD(SSD.INFO, SetIndex)
	re<-SKAT_emmaX(Z, obj, ...)
	return(re)
}


#
# Only SKAT_Null_Model obj can be used
#
SKAT_emmaX.SSD.All = function(SSD.INFO, obj, ...){
	
	N.Set<-SSD.INFO$nSets
	OUT.Pvalue<-rep(NA,N.Set)
	OUT.Marker<-rep(NA,N.Set)
	OUT.Marker.Test<-rep(NA,N.Set)
	OUT.Error<-rep(-1,N.Set)
	OUT.Pvalue.Resampling<-NULL

	Is.Resampling = FALSE
	n.Resampling = 0
	
	for(i in 1:N.Set){
		Is.Error<-TRUE
		try1<-try(SKAT:::Get_Genotypes_SSD(SSD.INFO, i),silent = TRUE)
		if(class(try1) != "try-error"){
			Z<-try1
			Is.Error<-FALSE
			
			
		} else {
			err.msg<-geterrmessage()
			msg<-sprintf("Error to get genotypes of %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
			warning(msg,call.=FALSE)
		}
	
		if(!Is.Error){
			Is.Error<-TRUE
			try2<-try(SKAT_emmaX(Z,obj, ...),silent = TRUE)
			
			if(class(try2) != "try-error"){
				re<-try2
				Is.Error<-FALSE
			} else {

				err.msg<-geterrmessage()
				msg<-sprintf("Error to run SKAT for %s: %s",SSD.INFO$SetInfo$SetID[i], err.msg)
				warning(msg,call.=FALSE)
			}
		}
		
		if(!Is.Error){

			OUT.Pvalue[i]<-re$p.value
			OUT.Marker[i]<-re$param$n.marker
			OUT.Marker.Test[i]<-re$param$n.marker.test
		}
	}

	
	out.tbl<-data.frame(SetID=SSD.INFO$SetInfo$SetID, P.value=OUT.Pvalue, N.Marker.All=OUT.Marker, N.Marker.Test=OUT.Marker.Test)
	re<-list(results=out.tbl)
	class(re)<-"SKAT_SSD_ALL"

	return(re)	
}

load_SKAT_NULL_emmaX<- function(formula, remlf, eigf, data = NULL) {
  #formula, data=NULL, K=NULL, Kin.File=NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
  eig <- .Call("readEigenWithIDs",eigf)
  remls <- read.table(remlf)[,2]
  delta <- remls[1]
  va <- remls[2]
  ve <- remls[3]
  maxLL <- remls[4]

  obj<-model.frame(formula,na.action = na.pass,data)

  n<-dim(obj)[1]
  id_include<-SKAT:::SKAT_Null_Model_Get_Includes(obj,obj)
  X<-model.matrix(formula,data=data)
  y<-model.response(obj)
  q <- ncol(X)

  P <- t(eig$evecs)
  res <- P %*% y

  re<-list( LL=maxLL, va=va, ve=ve, P=P, res=res, id_include=id_include )
  class(re)<-"SKAT_NULL_Model_EMMAX"
  return (re)
}

## SKAT
## binary vs quantitiative
## weighted vs flat
## optimal vs original
group.mmskat <- function() {
  ## flip ref/alt to minor/major
  flip <- which(AC>NS)
  genos[flip,] <- 2-genos[flip,]
  m <- nrow(genos)
  cname <- c("QSTAT");

  genos.nona <- genos
  genos.mean <- rowMeans(genos,na.rm=TRUE)  
  genos.mean.matrix <- matrix(genos.mean,nrow(genos),ncol(genos))
  genos.nona[is.na(genos)] <- genos.mean.matrix[is.na(genos)]  ## fill in missing genotypes
  genos.maf <- genos.mean/2
  if ( qr(genos.nona)$rank < 1 ) {  ## Change suggested by Jason
    m <- 0
  }  
  if ( m > 0 ) {
    eig <- .Call("readEigenWithIDs",eigf)
    remls <- as.matrix(read.table(remlf)[,2])
    delta <- remls[1]
    emxT <- matrix( 1/sqrt(eig$evals+delta), ncol(eig$evecs), nrow(eig$evecs) ) * t(eig$evecs)  ## not matrix multiplication
    tpheno <- emxT %*% pheno 
    tgenos <- emxT %*% t(genos.nona)

    obj <- SKAT_Null_Model(tpheno~1,out_type="C",n.Resampling=0,type.Resampling="bootstrap",Adjustment=skatAdjust)

    if ( skatOptimal ) {
      if ( skatFlat ) {
        if ( skatAdjust ) {
          r <- SKAT(tgenos,obj,kernel="linear",method="optimal.adj",impute.method="fixed", emmax_maf = genos.maf)
        }
        else {
          r <- SKAT(tgenos,obj,kernel="linear",method="optimal",impute.method="fixed", emmax_maf = genos.maf)
        }
      }
      else {
        if ( skatAdjust ) {
          r <- SKAT(tgenos,obj,kernel="linear.weighted",method="optimal.adj",weights.beta=betas,impute.method="fixed", emmax_maf = genos.maf)
        }
        else {
          r <- SKAT(tgenos,obj,kernel="linear.weighted",method="optimal",weights.beta=betas,impute.method="fixed", emmax_maf = genos.maf)
        }
      }

      p <- r$p.value
      if ( m == 1 ) {
        rho <- -998 # single variant
      }
      else {
        rho <- r$param$rho_est
        if ( length(rho) > 1 ) {
          rho <- -999  ## multiple rho temporary solution
        }
        if ( is.null(rho) ) {
          rho <- -997  ## rho is NULL temporary solution
        }
        
      }
      return(list(p=p,add=c(rho),cname=cname))
    }
    else {
      if ( skatFlat ) {
        r <- SKAT(tgenos,obj,kernel="linear",method="davies",impute.method="fixed", emmax_maf = genos.maf)
      }
      else {
        r <- SKAT(tgenos,obj,kernel="linear.weighted",method="davies",weights.beta=betas,impute.method="fixed", emmax_maf = genos.maf)
        #stop("foo")                
      }
      Q <- r$Q
      p <- ifelse(r$p.value==0,r$param$liu_pval,r$p.value)
      rm(r)
      gc()
      return(list(p=p,add=c(Q),cname=cname))
    }    
    #r <- SKAT(tgenos,obj,kernel="linear",method="davies",impute.method="fixed")
    #Q <- r$Q
    #p <- ifelse(r$p.value==0,p$param$liu <- pval,r$p.value)
    #rm(r)
    #gc()
    #return(list(p=p,add=c(Q),cname=cname))
  }
  else {
    return(list(p=NA,add=rep(NA,length(cname)),cname=cname))
  }
}


