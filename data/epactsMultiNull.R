args <- commandArgs(trailingOnly=TRUE)
nullf <- args[2]  ## output name
phenof <- args[3]  ## phenotype file
covf <- args[4]    ## covariate file or null model index file

pnames <- scan(phenof,what=character(),nlines=1)[-1]
phenos <- as.matrix(read.table(phenof)[,-1])

if ( covf == "NULL" ) {
  cov <- NULL;
} else {
  cov <- as.matrix(read.table(covf)[,-1])
}



ScoreTest_wSaddleApprox_Get_X1 = function(X1)
{
	q1<-ncol(X1)
	if(q1>=2)
	{
		if(sum(abs(X1[,1]-X1[,2]))==0)
		{
			X1=X1[,-2]
			q1<-q1-1
		}
	}
	qr1<-qr(X1)
	if(qr1$rank < q1){
		
		X1.svd<-svd(X1)
		X1 = X1.svd$u[,1:qr1$rank]
	} 

	return(X1)
}


ScoreTest_wSaddleApprox_NULL_Model <- function(formula, data=NULL)
{
	X1<-model.matrix(formula,data=data)
	X1<-ScoreTest_wSaddleApprox_Get_X1(X1)
	
	glmfit= glm(formula, data=data, family = "binomial")
  	mu    = glmfit$fitted.values
	V = mu*(1-mu)
  	res = glmfit$y- mu
	n1<-length(res)
	
	XV = t(X1 * V)
	XVX_inv= solve(t(X1)%*%(X1 * V))
	XXVX_inv= X1 %*% XVX_inv   
	
	re<-list(y=glmfit$y, cov=X1, mu=mu, res=res, V=V, X1=X1, XV=XV, XXVX_inv =XXVX_inv)
	class(re)<-"SA_NULL"
	return(re)	
}

##################################################################
## MAIN FUNCTION:  null.b.sna2
##################################################################

phenos <- phenos - min(phenos,na.rm=T)

    k <- ncol(cov)

g<-ncol(phenos)

	for(gind in 1:g)
	{
	nm<-which(is.na(phenos[,gind])==FALSE)
	pheno<-phenos[nm,gind]
	cov1<-cov[nm,]

	obj.null<-ScoreTest_wSaddleApprox_NULL_Model(pheno ~as.matrix(cov1))
        save.image(paste(nullf,".",pnames[gind],".null.RData",sep=""))

   	 }


