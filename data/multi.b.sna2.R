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
	
	re<-list(y=glmfit$y, mu=mu, res=res, V=V, X1=X1, XV=XV, XXVX_inv =XXVX_inv)
	class(re)<-"SA_NULL"
	return(re)	
}

Korg<-function(t, mu, g)
{
	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp<-log(1 - mu + mu * exp(g* t1))
		out[i]<-sum(temp)
	}
	return(out)
}


K1_adj<-function(t, mu, g, q)
{
	n.t<-length(t)	
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp1<-(1 - mu)* exp(-g * t1) + mu
		temp2<-mu *g
		out[i]<-sum(temp2/temp1)-q
	}
	return(out)
}

K2<-function(t, mu, g)
{
	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp1<-((1 - mu)* exp(-g * t1) + mu)^2
		temp2<-(1-mu) * mu * g^2 * exp(-g*t1)
		out[i]<-sum(temp2/temp1,na.rm=TRUE)
	}
	return(out)
}

getroot_K1<-function(init,mu,g,q,m1,tol=.Machine$double.eps^0.25,maxiter=1000)
{
	g.pos<-sum(g[which(g>0)])
	g.neg<- sum(g[which(g<0)])
	if(q>=g.pos || q<=g.neg)
	{
		return(list(root=Inf,n.iter=0,Is.converge=TRUE))
	} else {
		t<-init
		K1_eval<-K1_adj(t,mu,g,q)
		prevJump<- Inf
		rep<-1
		repeat
		{
			K2_eval<-K2(t,mu,g)
			tnew<-t-K1_eval/K2_eval
			if(is.na(tnew))
			{
				conv=FALSE
				break
			}
			if(abs(tnew-t)<tol)
			{
				conv<-TRUE
				break
			}
			if(rep==maxiter)
			{
				conv<-FALSE
				break
			}

			newK1<-K1_adj(tnew,mu,g,q)
			if(sign(K1_eval)!=sign(newK1))
			{
				if(abs(tnew-t)>prevJump-tol)
				{
					tnew<-t+sign(newK1-K1_eval)*prevJump/2
					newK1<-K1_adj(tnew,mu,g,q)
					prevJump<-prevJump/2
				} else {
					prevJump<-abs(tnew-t)
				}
			}

			rep<-rep+1
			t<-tnew
			K1_eval<-newK1
		}
		return(list(root=t,n.iter=rep,Is.converge=conv))
	}
}

Get_Saddle_Prob<-function(zeta, mu, g, q) 
{
	k1<-Korg(zeta, mu, g)
	k2<-K2(zeta, mu, g)
	
	if(is.finite(k1) && is.finite(k2))
	{
	temp1<-zeta * q - k1

	
	w<-sign(zeta) * (2 *temp1)^{1/2}
	v<- zeta * (k2)^{1/2}
			
	Z.test<-w + 1/w * log(v/w)	
	
	if(Z.test > 0){
		pval<-pnorm(Z.test, lower.tail = FALSE)
	} else {
		pval<-pnorm(Z.test, lower.tail = TRUE)
	}	
	} else {
	pval<-0
	}
	
	return(pval)
}
	
Saddle_Prob<-function(q, mu, g, Cutoff=2,alpha)
{
	m1<-sum(mu * g)
	var1<-sum(mu * (1-mu) * g^2)
	p1=NULL
	p2=NULL

	#
	qinv = -sign(q-m1) * abs(q-m1) + m1

	# Noadj
	pval.noadj<-pchisq((q - m1)^2/var1, lower.tail = FALSE, df=1)
	Is.converge=TRUE

 	if(Cutoff=="BE"){
		rho<-sum(((abs(g))^3)*mu*(1-mu)*(mu^2+(1-mu)^2))
		B<-0.56*rho*var1^(-3/2)
		p<-B+alpha/2
		Cutoff=ifelse(p>=0.496,0.01,qnorm(p,lower.tail=F))
	} else if(Cutoff < 10^-1){
		Cutoff=10^-1
	} 			
	#

	if(abs(q - m1)/sqrt(var1) < Cutoff || abs(-q - m1)/sqrt(var1) < Cutoff){

		pval=pval.noadj
		
	} else {
		out.uni1<-getroot_K1(0, mu=mu, g=g, q=q)
		out.uni2<-getroot_K1(0, mu=mu, g=g, q=qinv)
		if(out.uni1$Is.converge==TRUE && out.uni2$Is.converge==TRUE)
		{
			p1<-Get_Saddle_Prob(out.uni1$root, mu, g, q)
			p2<-Get_Saddle_Prob(out.uni2$root, mu, g, qinv)
			pval = p1+p2
			Is.converge=TRUE
		} else {
 			print("Error_Converge")
			pval<-pval.noadj
			Is.converge=FALSE	
		}				
	}
	
	return(list(p.value=pval, p.value.NA=pval.noadj, Is.converge=Is.converge, p1=p1, p2=p2))
}

TestSPA<-function(G, obj.null, Cutoff=2,alpha)
{
	if(class(obj.null) != "SA_NULL"){
		stop("obj.null should be a returned object from ScoreTest_wSaddleApprox_NULL_Model")
	}
	
	y = obj.null$y
	mu = obj.null$mu
	res = obj.null$res
	
	n.g<-sum(G)
	if(n.g/(2*length(G))>0.5)
	{
		G<-2-G
		n.g<-sum(G)
	}
	G1<-G  -  obj.null$XXVX_inv %*%  (obj.null$XV %*% G)
	q<-sum(G1 * y) /sqrt(n.g)
	out<-Saddle_Prob(q, mu=mu, g=G1/sqrt(n.g), Cutoff=Cutoff,alpha=alpha)

	return(out)
}

Korg_fast<-function(t, mu, g,gNA,gNB,muNA,muNB,NAmu,NAsigma)
{
	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp<-log(1 - muNB + muNB * exp(gNB* t1))
		out[i]<-sum(temp)+NAmu*t1+0.5*NAsigma*t1^2
	}
	return(out)
}


K1_adj_fast<-function(t, mu, g, q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
{
	n.t<-length(t)	
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp1<-(1 - muNB)* exp(-gNB * t1) + muNB
		temp2<-muNB *gNB
		temp3<-NAmu+NAsigma*t1
		out[i]<-sum(temp2/temp1)+temp3-q
	}
	return(out)
}

K2_fast<-function(t, mu, g,gNA,gNB,muNA,muNB,NAmu,NAsigma)
{
	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp1<-((1 - muNB)* exp(-gNB * t1) + muNB)^2
		temp2<-(1-muNB) * muNB * gNB^2 * exp(-gNB*t1)
		out[i]<-sum(temp2/temp1,na.rm=TRUE)+NAsigma
	}
	return(out)
}

getroot_K1_fast<-function(init,mu,g,q,m1,gNA,gNB,muNA,muNB,NAmu,NAsigma,tol=.Machine$double.eps^0.25,maxiter=1000)
{
	g.pos<-sum(g[which(g>0)])
	g.neg<- sum(g[which(g<0)])
	if(q>=g.pos || q<=g.neg)
	{
		return(list(root=Inf,n.iter=0,Is.converge=TRUE))
	} else {
		t<-init
		K1_eval<-K1_adj_fast(t,mu,g,q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
		prevJump<- Inf
		rep<-1
		repeat
		{
			K2_eval<-K2_fast(t,mu,g,gNA,gNB,muNA,muNB,NAmu,NAsigma)
			tnew<-t-K1_eval/K2_eval
			if(is.na(tnew))
			{
				conv=FALSE
				break
			}
			if(abs(tnew-t)<tol)
			{
				conv<-TRUE
				break
			}
			if(rep==maxiter)
			{
				conv<-FALSE
				break
			}

			newK1<-K1_adj_fast(tnew,mu,g,q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
			if(sign(K1_eval)!=sign(newK1))
			{
				if(abs(tnew-t)>prevJump-tol)
				{
					tnew<-t+sign(newK1-K1_eval)*prevJump/2
					newK1<-K1_adj_fast(tnew,mu,g,q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
					prevJump<-prevJump/2
				} else {
					prevJump<-abs(tnew-t)
				}
			}

			rep<-rep+1
			t<-tnew
			K1_eval<-newK1
		}
		return(list(root=t,n.iter=rep,Is.converge=conv))
	}
}

Get_Saddle_Prob_fast<-function(zeta, mu, g, q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
{
	k1<-Korg_fast(zeta, mu, g,gNA,gNB,muNA,muNB,NAmu,NAsigma)
	k2<-K2_fast(zeta, mu, g,gNA,gNB,muNA,muNB,NAmu,NAsigma)
	
	if(is.finite(k1) && is.finite(k2))
	{
	temp1<-zeta * q - k1

	
	w<-sign(zeta) * (2 *temp1)^{1/2}
	v<- zeta * (k2)^{1/2}
			
	Z.test<-w + 1/w * log(v/w)	
	
	if(Z.test > 0){
		pval<-pnorm(Z.test, lower.tail = FALSE)
	} else {
		pval<-pnorm(Z.test, lower.tail = TRUE)
	}	
	} else {
	pval<-0
	}
	
	return(pval)
}

Saddle_Prob_fast<-function(q, g,mu,gNA,gNB,muNA,muNB,Cutoff=2,alpha)
{
	m1<-sum(mu * g)
	var1<-sum(mu * (1-mu) * g^2)
	p1=NULL
	p2=NULL

	#
	qinv = -sign(q-m1) * abs(q-m1) + m1

	# Noadj
	pval.noadj<-pchisq((q - m1)^2/var1, lower.tail = FALSE, df=1)
	Is.converge=TRUE

	if(Cutoff=="BE"){
		rho<-sum(((abs(g))^3)*mu*(1-mu)*(mu^2+(1-mu)^2))
		B<-0.56*rho*var1^(-3/2)
		p<-B+alpha/2
		Cutoff=ifelse(p>=0.496,0.01,qnorm(p,lower.tail=F))
	} else if(Cutoff < 10^-1){
		Cutoff=10^-1
	} 
			
	#
	
	if(abs(q - m1)/sqrt(var1) < Cutoff || abs(-q - m1)/sqrt(var1) < Cutoff){

		pval=pval.noadj
		
	} else {
		NAmu= m1-sum(gNB*muNB)
		NAsigma=var1-sum(muNB*(1-muNB)*gNB^2)
		out.uni1<-getroot_K1_fast(0, mu=mu, g=g, q=q,gNA=gNA,gNB=gNB,muNA=muNA,muNB=muNB,NAmu=NAmu,NAsigma=NAsigma)
		out.uni2<-getroot_K1_fast(0, mu=mu, g=g, q=qinv,gNA=gNA,gNB=gNB,muNA=muNA,muNB=muNB,NAmu=NAmu,NAsigma=NAsigma)
		if(out.uni1$Is.converge==TRUE && out.uni2$Is.converge==TRUE)
		{
			p1<-Get_Saddle_Prob_fast(out.uni1$root, mu, g, q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
			p2<-Get_Saddle_Prob_fast(out.uni2$root, mu, g, qinv,gNA,gNB,muNA,muNB,NAmu,NAsigma)
			pval = p1+p2
			Is.converge=TRUE
		} else {
 			print("Error_Converge")
			pval<-pval.noadj
			Is.converge=FALSE	
		}		
		
	}
	
	return(list(p.value=pval, p.value.NA=pval.noadj, Is.converge=Is.converge, p1=p1, p2=p2))
}

TestSPAfast<-function(G, obj.null, Cutoff=2,alpha)
{
	if(class(obj.null) != "SA_NULL"){
		stop("obj.null should be a returned object from ScoreTest_wSaddleApprox_NULL_Model")
	}
	
	y = obj.null$y
	mu = obj.null$mu
	res = obj.null$res
	
	n.g<-sum(G)
	if(n.g/(2*length(G))>0.5)
	{
		G<-2-G
		n.g<-sum(G)
	}
	NAset<-which(G==0)
	G1<-G  -  obj.null$XXVX_inv %*%  (obj.null$XV %*% G)
	q<-sum(G1 * y) /sqrt(n.g)
	g=G1/sqrt(n.g)
	if(length(NAset)/length(G)<0.5)
	{
		out<-Saddle_Prob(q, mu=mu, g=G1/sqrt(n.g), Cutoff=Cutoff,alpha=alpha)

	} else {

	out<-Saddle_Prob_fast(q,g=g,mu=mu,gNA=g[NAset],gNB=g[-NAset],
muNA=mu[NAset],muNB=mu[-NAset],Cutoff=Cutoff,alpha=alpha)
	}
	return(out)
}


##################################################################
## MAIN FUNCTION:  multi.b.sna2
##################################################################
multi.b.sna2 <- function() {
	#write.table(pheno,"y",row.names=F,col.names=F)
	#write.table(genos,"genos",row.names=F,col.names=F)
	#write.table(cov,"cov",row.names=F,col.names=F)
phenos <- phenos - min(phenos,na.rm=T)

    n <- ncol(genos)
    k <- ncol(cov)
	genos<-as.matrix(genos)
	if(ncol(genos)==1)
	{
		m<-1
		genos<-t(genos)
	} else {
		m <- nrow(genos)
	}
g<-ncol(phenos)

    ## resolve missing genotype by mean imputation
    ina <- is.na(genos)
    if ( length(vids) > 0 ) {
        genos[ina] <- matrix(AC[vids]/NS[vids],nrow(genos),ncol(genos))[ina]
    }
    
    p <- matrix(NA,m,g) # store p-values for m markers
    add <- matrix(NA,m,g) # extra columns:  Beta, SE, Chisq
    cname <- paste(pnames,"P.NA",sep=".") #,"N.CASE","N.CTRL","AF.CASE","AF.CTRL")

	for(gind in 1:g)
	{
    if ( m > 0 ) { ## If there is at least one marker to test
	nm<-which(is.na(phenos[,gind])==FALSE)
	pheno<-phenos[nm,gind]
	cov1<-cov[nm,]
	geno<-as.matrix(genos[,nm])
	if(ncol(geno)==1)	geno<-t(geno)

	obj.null<-ScoreTest_wSaddleApprox_NULL_Model(pheno ~as.matrix(cov1))
        
        for (i in 1:m) {
            re <- TestSPAfast(as.vector(geno[i,,drop=FALSE]), obj.null, Cutoff=2)
            p[i,gind] <- re$p.value
            
            add[i,gind] <- re$p.value.NA
        }
    }
	}


    return(list(p=p, add=add, cname=cname))
}


