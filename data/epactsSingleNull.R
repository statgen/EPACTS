args <- commandArgs(trailingOnly=TRUE)
nullf <- args[2]  ## output name
phenof <- args[3]  ## phenotype file
covf <- args[4]    ## covariate file or null model index file

pheno <- as.matrix(read.table(phenof)[,-1])

if ( covf == "NULL" ) {
  cov <- NULL;
} else {
  cov <- as.matrix(read.table(covf)[,-1])
}


fast.logistf.fit <- function (x, y, weight = NULL, offset = NULL, firth = TRUE, col.fit = NULL, 
    init = NULL, control) {
    n <- nrow(x)
    k <- ncol(x)
    if (is.null(init)) 
        init = rep(0, k)
    if (is.null(col.fit)) 
        col.fit = 1:k
    if (is.null(offset)) 
        offset = rep(0, n)
    if (is.null(weight)) 
        weight = rep(1, n)
    if (col.fit[1] == 0) 
        maxit <- 0
    if (missing(control)) 
        control <- fast.logistf.control()
    maxit <- control$maxit
    maxstep <- control$maxstep
    maxhs <- control$maxhs
    lconv <- control$lconv
    gconv <- control$gconv
    xconv <- control$xconv
    beta <- init
    iter <- 0
    pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))
    evals <- 1
    repeat {
        beta.old <- beta
        XW2 <- t(x * (weight * pi * (1-pi))^0.5)
        myQR <- qr(t(XW2))
        Q <- qr.Q(myQR)
        h <- (Q*Q) %*% rep(1, ncol(Q))        
        if (firth) 
            U.star <- crossprod(x, weight * (y - pi) + h * (0.5 - pi))
        else U.star <- crossprod(x, weight * (y - pi))
        XX.covs <- matrix(0, k, k)
        if (col.fit[1] != 0) {
            XX.XW2 <- t(x[, col.fit, drop=FALSE] * (weight * pi * (1-pi))^0.5)
            XX.Fisher <- crossprod(t(XX.XW2))
            XX.covs[col.fit, col.fit] <- fast.invFisher(XX.Fisher)   ###### HERE IS THE PROBLEM!!!
        }
        if(all(is.na(XX.covs)) == T) {
            break
        }  
        delta <- as.vector(XX.covs %*% U.star)
        delta[is.na(delta)] <- 0
        mx <- max(abs(delta))/maxstep
        if (mx > 1) 
            delta <- delta/mx
        evals <- evals + 1
        if (maxit > 0) {
            iter <- iter + 1
            beta <- beta + delta
                pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))


        }
        if (iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star[col.fit]) < 
            gconv)))) 
            break
    }
    # Error catching (if chol(x) not positive definite)
    if(all(is.na(XX.covs))==T) {
        var <- XX.covs
        list(beta = NA, var = var, pi = NA, hat.diag = NA, 
  iter = NA, evals = NA, conv = c(NA, 
            NA, NA))
    } else {
        var <- XX.covs
        list(beta = beta, var = var, pi = pi, hat.diag = h, 
            iter = iter, evals = evals, conv = c(max(abs(U.star)),
		 max(abs(delta))))
    }
}


fast.logistf.control <- function (maxit = 50, maxhs = 15, maxstep = 15, lconv = 1e-05, 
    gconv = 1e-05, xconv = 1e-05) 
{
    list(maxit = maxit, maxhs = maxhs, maxstep = maxstep, lconv = lconv, 
        gconv = gconv, xconv = xconv)
}

fast.logDet <- function (x) {
    my.chol <- tryCatch(chol(x),error=function(e) {NA})
    if (all(is.na(my.chol))==T) {
        return(NA)
    } else {
        return (2 * sum(log(diag(my.chol))))
    }
}   

fast.invFisher <- function(x) {
  my.chol <- tryCatch(chol(x),error=function(e) {NA})
    if (all(is.na(my.chol))==T) {
        return(NA)
    } else {
        return (chol2inv(my.chol))
    }
  #ifelse(is.na(my.chol), NA, chol2inv(my.chol))
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
	convflag=0
	if(glmfit$converged)
	{
		mu = glmfit$fitted.values
		if(mean(mu)/mean(glmfit$y)>0.001 & (1-mean(mu))/(1-mean(glmfit$y))>0.001)	convflag<-1	#Check that the null model converged properly with glm
	}
	if(convflag==0)
	{
		firthfit=fast.logistf.fit(x=X1,y=glmfit$y)
		eta<-X1%*%firthfit$beta
  		mu    = as.vector(exp(eta)/(1+exp(eta)))
	}

	V = mu*(1-mu)
  	res = glmfit$y- mu
	n1<-length(res)
	
	XV = t(X1 * V)
	XVX_inv= solve(t(X1)%*%(X1 * V))
	XXVX_inv= X1 %*% XVX_inv   
	
	re<-list(y=glmfit$y, mu=mu, res=res, V=V, X1=X1, XV=XV, XXVX_inv =XXVX_inv)
	class(re)<-"SA_NULL"
	return(re)	
}

##################################################################
## MAIN FUNCTION:  null.b.sna2
##################################################################

pheno <- pheno - min(pheno,na.rm=T)
    k <- ncol(cov)

	obj.null<-ScoreTest_wSaddleApprox_NULL_Model(pheno ~as.matrix(cov))
        save.image(paste(nullf,".null.RData",sep=""))




