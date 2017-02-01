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

## single.b.firth :
## Use Firth bias-reduced logistic regression to perform association
## By: Clement Ma
## 
## Adapted from 'logistf' R package (v1.10)
## By:  Ploner M, Dunkler D, Southworth H, Heinze G
## TRAITS  : BINARY
## RETURNS : PVALUE, BETA, SEBETA, CHISQ
## MISSING VALUES : IGNORED

Korg<-function(t, mu, g){

	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp<-log(1 - mu + mu * exp(g* t1))
		out[i]<-sum(temp)
	}
	return(out)
}


K1_adj<-function(t, mu, g, q){

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


K2<-function(t, mu, g){

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

Get_Spline_Fun<-function(mu, g){
	node<-c(-5000,-2000,-1000, -100, -50.5, -38.125, -25.75, -13.375, -1, 0,
	1, 50.5, 62.875, 75.25, 100, 1000,2000,5000)
	y1<-Korg(node ,mu ,g)
	y1deriv<-K1_adj(node ,mu ,g,0)
	sfun<-splinefun(node , y1 ,method='natural')

}

Get_Spline_FunH<-function(mu, g){
	
	node<-c(-5000,-2000,-1000, -100, -1, 0, 1, 100, 1000,2000,5000)
	flag=0
	rep<-0
	repeat{
		y1<-Korg(node ,mu ,g)
		y1deriv<-K1_adj(node ,mu ,g,0)
		sfun<-splinefunH(node , y1 ,m=y1deriv)
		xval<-NULL
		for(i in 1:(length(node)-1))
			xval<-c(xval,seq(node[i],node[i+1],length.out=100))
		curve<-sfun(xval,deriv=2)
		if(min(curve)<0)
		{
			ls<-unique(findInterval(xval[which(curve<0)],node,all.inside=T))
			nodeadd<-NULL
			for(i in 1:length(ls))
			{
				b<-(node[ls[i]]+node[ls[i]+1])/2
				nodeadd<-c(nodeadd,b)	
			}
			nodeadd<-unique(nodeadd)
			node<-c(node,nodeadd)
			node<-node[order(node)]
		} else {
			flag=1
		}
		if(flag==1 || rep>10) break
		rep<-rep+1
	}


	return(sfun)
}

Get_Spline_Fun_H2<-function(mu, g){
	
	node<-c(-5000,-2000,-1000, -100, -1, 0, 1, 100, 1000,2000,5000)
	y1<-Korg(node ,mu ,g)
	y1deriv<-K1_adj(node ,mu ,g,0)
	y1deriv2<-K2(node ,mu ,g)
	sfun<-splineH2(nodeval=cbind(node,y1,y1deriv,y1deriv2))
	return(sfun)
}

splineH2<-function(nodeval)
{
	sfun<-matrix(0,nrow(nodeval),7)
	sfun[,1]<-nodeval[,1]
	A<-matrix(c(1,1,0,0,0,0,0,1,1,1,0,0,0,1,0,2,2,2,0,1,0,3,0,6,0,1,0,
	4,0,12,0,1,0,5,0,20),6,6)
	for(i in 1:(nrow(nodeval)-1))
	{
		f<-as.vector(nodeval[i:(i+1),2:4])
		sfun[i,2:7]<-solve(A)%*%f
	}
	return(sfun)
}

Korg_S_H2<-function(t, sfun){

	i<-findInterval(t,sfun[,1],all.inside=T)
	t1<-(t-sfun[i,1])/(sfun[i+1,1]-sfun[i,1])
	x<-t1^(0:5)
	out<-sum(x*sfun[i,2:7])
	return(out)
}


K1_adj_S_H2<-function(t, q, sfun){

	i<-findInterval(t,sfun[,1],all.inside=T)
	t1<-(t-sfun[i,1])/(sfun[i+1,1]-sfun[i,1])
	x<-(t1^(0:4))*(1:5)
	out<-sum(x*sfun[i,3:7])-q
	return(out)
}


K2_S_H2<-function(t, sfun){

	i<-findInterval(t,sfun[,1],all.inside=T)
	t1<-(t-sfun[i,1])/(sfun[i+1,1]-sfun[i,1])
	x<-(t1^(0:3))*c(2,6,12,20)
	out<-sum(x*sfun[i,4:7])
	return(out)
}

Korg_S<-function(t, sfun){

	out<-sfun(t, deriv=0)
	return(out)
}


K1_adj_S<-function(t, q, sfun){

	out<-sfun(t, deriv=1) - q
	return(out)
}


K2_S<-function(t, sfun){

	out<-sfun(t, deriv=2) 
	return(out)
}


Get_Saddle_Prob<-function(zeta, mu, g, q) {
	
	# zeta<-out.uni2$root
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
	


Saddle_Prob<-function(q, mu, g, Cutoff=2){

	#mu =pi; g=g=G1/sqrt(n.g); Spline=FALSE
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
		alpha<-5*10^-5
		p<-B+alpha/2
		Cutoff=ifelse(p>=0.496,0.01,qnorm(p,lower.tail=F))
	} else if(Cutoff < 10^-1){
		Cutoff=10^-1
	} 
	
	#
	
	if(abs(q - m1)/sqrt(var1) < Cutoff || abs(-q - m1)/sqrt(var1) < Cutoff){
		#Ques: Why the 2nd condition? Isn't that redundant as q and m1 are both positives,
		#and in general can cause problem if q and m1 has different signs?

		pval=pval.noadj
		
	} else {
		out.uni1<-try(uniroot(K1_adj, c(-100,100), mu=mu, g=g, q=q), silent=TRUE)
		out.uni2<-try(uniroot(K1_adj, c(-100,100), mu=mu, g=g, q=qinv), silent=TRUE)
		
#Ques: Don't need the whole interval
		# try one-more with large interval
		if(class(out.uni1) == "try-error"){
			out.uni1<-try(uniroot(K1_adj, c(-1000,1000), mu=mu, g=g, q=q), silent=TRUE)
		}		
		
		if(class(out.uni1) != "try-error"){

			p1<-Get_Saddle_Prob(out.uni1$root, mu, g, q)
			
			if(class(out.uni2) != "try-error"){
				p2<-Get_Saddle_Prob(out.uni2$root, mu, g, qinv)
				
			} else {
				t<-seq(-1000,1000, length.out=20)
				out1<-K1_adj(t, mu, g, qinv)
				idx<-which(abs(out1) == min(abs(out1), na.rm=TRUE))
				try(p2<-Get_Saddle_Prob(t[idx], mu, g, qinv), silent=TRUE)
				if(class(p2) == "try-error"){
					p2  = p1
				}
				
			}
			
			pval = p1+p2
			
			#if(pval < 10^-3 && pval.noadj > 0.1){
			#	pval<-pval.noadj
			#	Is.converge=FALSE				
			#}
		
		} else {
			pval<-pval.noadj
			Is.converge=FALSE
		}
		
		
	} 
	
	return(list(p.value=pval, p.value.NA=pval.noadj, Is.converge=Is.converge, p1=p1, p2=p2))
}

Saddle_Prob_NR<-function(q, mu, g, Cutoff=2){

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
		alpha<-5*10^-5
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
		#if(pval < 10^-3 && pval.noadj > 0.1){
 		#	print("Error_Mismatch")
		#	pval<-pval.noadj
		#	Is.converge=FALSE				
		#}
		
	}
	
	return(list(p.value=pval, p.value.NA=pval.noadj, Is.converge=Is.converge, p1=p1, p2=p2))
}



Saddle_Prob_wSpline<-function(q, mu, g, Spline=FALSE){

	#mu =pi; g=g=G1/sqrt(n.g); Spline=FALSE
	m1<-sum(mu * g)
	var1<-sum(mu * (1-mu) * g^2)

	# Noadj
	pval.noadj<-pchisq((q - m1)^2/var1, lower.tail = FALSE, df=1)
	pval.noadj<-pnorm((q - m1)/sqrt(var1), lower.tail = FALSE)
	Is.converge=TRUE
		
	# incase of singleton
	idx<-which(g> 0)		#Ques: Why these steps?
	if(length(idx)==1){
		if(q==0){
			pval=1 - (1-mu[idx])/2
		} else {
			pval=mu[idx]/2
		}
	} else if(q ==0){
		
		temp<-sum(log(1-mu[idx]))
		pval<-1 - exp(temp)/2
		
	} else if (q==sum(g)){
		
		temp<-sum(log(mu[idx]))
		pval<-exp(temp)/2
		
	} else if(abs(q - m1)/sqrt(var1) < 10^{-5}){
		pval=0.5
	} else if(Spline==FALSE){
		
		out<-try(uniroot(K1_adj, c(-100,100), mu=mu, g=g, q=q), silent=TRUE)
		if(class(out) != "try-error"){
			zeta<-out$root
	
			k1<-Korg(zeta, mu, g)
			k2<-K2(zeta, mu, g)
	
			temp1<-zeta * q - k1

	
			w<-sign(zeta) * (2 *temp1)^{1/2}
			v<- zeta * (k2)^{1/2}
			
			Z.test<-w + 1/w * log(v/w)
			#pval<-pchisq(Z.test^2, lower.tail = FALSE, df=1)
			pval<-pnorm(Z.test, lower.tail = FALSE)
			
			if(pval < 10^-2 && pval.noadj > 0.1){
				pval<-pval.noadj
				Is.converge=FALSE				
			}
		
		} else {
			pval<-pval.noadj
			Is.converge=FALSE
		}
		
		
	} else {
		n<-length(g)
		sfun<-Get_Spline_Fun(mu, g, n)
		
		out<-uniroot(K1_adj_S, c(-3000,3000), q=q, sfun=sfun)
		zeta<-out$root
	
		k1<-Korg_S(zeta, sfun)
		k2<-K2_S(zeta, sfun)
	
		temp1<-zeta * q - k1

	
		w<-sign(zeta) * (2 *temp1)^{1/2}
		
		if(k2 <= 0){
			v<-w
		} else {
			v<- zeta * (k2)^{1/2}
		}
		pval<-pnorm(w + 1/w * log(v/w), lower.tail = FALSE)
	
	
	}
	q.transe<-q - m1
	if(!is.na(pval)){
		q.transe<-qnorm(pval, lower.tail=FALSE) * sqrt(var1) 
	} #Ques: Why this overwriting and why is this information needed

	
	return(list(pval=pval, pval.noadj=pval.noadj, q.transe=q.transe, Is.converge=Is.converge, w=w, v=v))
}


#
# Input
#	G: genotype
#	obj.null: returned object from GLM function with the null model 
#	Cutoff: a cutoff value to use either saddle point approx or normal approx
#
# Output
#	p.value: p-value from the saddle point approx. If the mean centered test statistic is within Cutoff *SD, p.value will be calculated using the normal approx.	
#	p.value.NA: p-value from the normal approx
#	Is.converge: covergence status of saddle point approx. If Is.converge=FALSE, p.value  will be calculated using the normal approx.
#	p1 and p2: internal use only (please ignore them)
#



#####################################################
# Following 3 functions are added

ScoreTest_wSaddleApprox_Get_X1 = function(X1){
	
	q1<-ncol(X1)
	if(q1>=3)
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
		X1 = X1.svd$u	#Ques: Why?
	} 

	return(X1)

}


ScoreTest_wSaddleApprox_NULL_Model = function(formula, data=NULL ){


	X1<-model.matrix(formula,data=data)
	X1<-ScoreTest_wSaddleApprox_Get_X1(X1)
	
	glmfit= glm(formula, data=data, family = "binomial")
  	mu    = glmfit$fitted.values
  	#Ques: Why not use Firth bias correction here? Not a big deal though, and not time consuming
	V = mu*(1-mu)
  	res = glmfit$y- mu
	n1<-length(res)
	
	#
	XV = t(X1 * V)
	XVX_inv= solve(t(X1)%*%(X1 * V))
	XXVX_inv= X1 %*% XVX_inv   
	
	re<-list(y=glmfit$y, mu=mu, res=res, V=V, X1=X1, XV=XV, XXVX_inv =XXVX_inv)
	class(re)<-"SA_NULL"
	return(re)
	
	
}


#
#	G: genotype vector
#	obj.SA.null: object from ScoreTest_wSaddleApprox_NULL_Model
#

ScoreTest_wSaddleApprox_New<-function(G, obj.SA.null, Cutoff=2){

	if(class(obj.SA.null) != "SA_NULL"){
		stop("obj.SA.null should be a returned object from ScoreTest_wSaddleApprox_NULL_Model")
	}
	
	y = obj.SA.null$y
	mu = obj.SA.null$mu
	res = obj.SA.null$res
	
	n.g<-sum(G)
	if(n.g/(2*length(G))>0.5)
	{
		G<-2-G
		n.g<-sum(G)
	}
	G1<-G  -  obj.SA.null$XXVX_inv %*%  (obj.SA.null$XV %*% G)
	#Ques: Reference? Doesn't removing X2 make G and y effectively indep?
	q<-sum(G1 * y) /sqrt(n.g)
	#q1<-sum(G * y) /sqrt(n.g)
	out<-Saddle_Prob(q, mu=mu, g=G1/sqrt(n.g), Cutoff=Cutoff)

	return(out)
}

ScoreTest_wSaddleApprox_NewNR<-function(G, obj.SA.null, Cutoff=2){

	if(class(obj.SA.null) != "SA_NULL"){
		stop("obj.SA.null should be a returned object from ScoreTest_wSaddleApprox_NULL_Model")
	}
	
	y = obj.SA.null$y
	mu = obj.SA.null$mu
	res = obj.SA.null$res
	
	n.g<-sum(G)
	if(n.g/(2*length(G))>0.5)
	{
		G<-2-G
		n.g<-sum(G)
	}
	G1<-G  -  obj.SA.null$XXVX_inv %*%  (obj.SA.null$XV %*% G)
	#Ques: Reference? Doesn't removing X2 make G and y effectively indep?
	q<-sum(G1 * y) /sqrt(n.g)
	#q1<-sum(G * y) /sqrt(n.g)
	out<-Saddle_Prob_NR(q, mu=mu, g=G1/sqrt(n.g), Cutoff=Cutoff)

	return(out)
}

ScoreTest_wSaddleApprox_OldNR<-function(G, obj.SA.null, Cutoff=2){

	if(class(obj.SA.null) != "SA_NULL"){
		stop("obj.SA.null should be a returned object from ScoreTest_wSaddleApprox_NULL_Model")
	}
	
	y = obj.SA.null$y
	mu = obj.SA.null$mu
	res = obj.SA.null$res
	
	n.g<-sum(G)
	if(n.g/(2*length(G))>0.5)
	{
		G<-2-G
		n.g<-sum(G)
	}
	NAset<-which(G==0)
	G1<-G  -  obj.SA.null$XXVX_inv %*%  (obj.SA.null$XV %*% G)
	q<-sum(G1 * y) /sqrt(n.g)
	g=G1/sqrt(n.g)
	if(length(NAset)/length(G)<0.5)
	{
		out<-Saddle_Prob_NR(q, mu=mu, g=G1/sqrt(n.g), Cutoff=Cutoff)

	} else {

	out<-Saddle_Prob_OldNR(q,g=g,mu=mu,gNA=g[NAset],gNB=g[-NAset],
muNA=mu[NAset],muNB=mu[-NAset],Cutoff=Cutoff)
	}
	return(out)
}

Saddle_Prob_OldNR<-function(q, g,mu,gNA,gNB,muNA,muNB,Cutoff=2){

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
		alpha<-5*10^-5
		p<-B+alpha/2
		Cutoff=ifelse(p>=0.496,0.01,qnorm(p,lower.tail=F))
	} else if(Cutoff < 10^-1){
		Cutoff=10^-1
	} 
			
	#
	
	if(abs(q - m1)/sqrt(var1) < Cutoff || abs(-q - m1)/sqrt(var1) < Cutoff){

		pval=pval.noadj
		
	} else {
		#NAmu= sum(gNA*muNA)
		#NAsigma=sum(muNA*(1-muNA)*gNA^2)
		NAmu= m1-sum(gNB*muNB)
		NAsigma=var1-sum(muNB*(1-muNB)*gNB^2)
		out.uni1<-getroot_K1_Old(0, mu=mu, g=g, q=q,gNA=gNA,gNB=gNB,muNA=muNA,muNB=muNB,NAmu=NAmu,NAsigma=NAsigma)
		out.uni2<-getroot_K1_Old(0, mu=mu, g=g, q=qinv,gNA=gNA,gNB=gNB,muNA=muNA,muNB=muNB,NAmu=NAmu,NAsigma=NAsigma)
		if(out.uni1$Is.converge==TRUE && out.uni2$Is.converge==TRUE)
		{
			p1<-Get_Saddle_Prob_Old(out.uni1$root, mu, g, q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
			p2<-Get_Saddle_Prob_Old(out.uni2$root, mu, g, qinv,gNA,gNB,muNA,muNB,NAmu,NAsigma)
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

getroot_K1_Old<-function(init,mu,g,q,m1,gNA,gNB,muNA,muNB,NAmu,NAsigma,tol=.Machine$double.eps^0.25,maxiter=1000)
{
	g.pos<-sum(g[which(g>0)])
	g.neg<- sum(g[which(g<0)])
	if(q>=g.pos || q<=g.neg)
	{
		return(list(root=Inf,n.iter=0,Is.converge=TRUE))
	} else {
		t<-init
		K1_eval<-K1_adj_Old(t,mu,g,q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
		prevJump<- Inf
		rep<-1
		repeat
		{
			K2_eval<-K2_Old(t,mu,g,gNA,gNB,muNA,muNB,NAmu,NAsigma)
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

			newK1<-K1_adj_Old(tnew,mu,g,q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
			if(sign(K1_eval)!=sign(newK1))
			{
				if(abs(tnew-t)>prevJump-tol)
				{
					tnew<-t+sign(newK1-K1_eval)*prevJump/2
					newK1<-K1_adj_Old(tnew,mu,g,q,gNA,gNB,muNA,muNB,NAmu,NAsigma)
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

Korg_Old<-function(t, mu, g,gNA,gNB,muNA,muNB,NAmu,NAsigma){

	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp<-log(1 - muNB + muNB * exp(gNB* t1))
		out[i]<-sum(temp)+NAmu*t1+0.5*NAsigma*t1^2
	}
	return(out)
}


K1_adj_Old<-function(t, mu, g, q,gNA,gNB,muNA,muNB,NAmu,NAsigma){

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

K2_Old<-function(t, mu, g,gNA,gNB,muNA,muNB,NAmu,NAsigma){

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

Get_Saddle_Prob_Old<-function(zeta, mu, g, q,gNA,gNB,muNA,muNB,NAmu,NAsigma) {
	
	# zeta<-out.uni2$root
	k1<-Korg_Old(zeta, mu, g,gNA,gNB,muNA,muNB,NAmu,NAsigma)
	k2<-K2_Old(zeta, mu, g,gNA,gNB,muNA,muNB,NAmu,NAsigma)
	
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

##################################################################
## MAIN FUNCTION:  single.b.sna0
##################################################################
single.b.sna0 <- function() {
    #control <- fast.logistf.control(maxit=25) # Max iterations to converge = 25

    n <- ncol(genos)
    k <- ncol(cov)
    m <- nrow(genos)

    ## resolve missing genotype by mean imputation
    ina <- is.na(genos)
    if ( length(vids) > 0 ) {
        genos[ina] <- matrix(AC[vids]/NS[vids],nrow(genos),ncol(genos))[ina]
    }
    
    p <- rep(NA,m) # store p-values for m markers
    add <- matrix(NA,m,2) # extra columns:  Beta, SE, Chisq
    cname <- c("PVAL.NA","CONVERGE") #,"N.CASE","N.CTRL","AF.CASE","AF.CTRL")

    if ( m > 0 ) { ## If there is at least one marker to test
	obj.SA.null<-ScoreTest_wSaddleApprox_NULL_Model(pheno ~as.matrix(cov))
        
        for (i in 1:m) {
            re <- ScoreTest_wSaddleApprox_OldNR(as.vector(genos[i,,drop=FALSE]), obj.SA.null, Cutoff=0)
            p[i] <- re$p.value
            
            add[i,1:2] <- c(re$p.value.NA, re$Is.converge)
        }
    }
	#write.table(pheno,"/net/wonderland/home/deyrnk/y",row.names=F,col.names=F)
	#write.table(genos,"/net/wonderland/home/deyrnk/genos",row.names=F,col.names=F)
	#write.table(cov,"/net/wonderland/home/deyrnk/cov",row.names=F,col.names=F)

    return(list(p=p, add=add, cname=cname))
}


