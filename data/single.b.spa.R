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

Get_Ks<-function(mu, g){


	mu1 = p *g
	mu2 = p * (1-p) *g^2
	mu3 = p * (1-p) * (1-2*p) *g^3
	mu4 = p * (1-p)* (3*p^2 - 3*p +1) *g^4
	
	k1 = mu1
	k2 = mu2
	k3 = mu3
	k4 = mu4 - 3* mu2^2

	return(cbind(k1,k2,k3,k4))

}


Korg_T4<-function(t, mu, g){

	ks<-Get_Ks(mu, g)
	k1 = ks[,1]
	k2 = ks[,2]
	k3 = ks[,3]
	k4 = ks[,4]
	
		
	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp<-k1* t1 + k2 * t1^2/2 + k3* t1^3/6 + k4 *t1^4 / 24
		out[i]<-sum(temp)
	}
	return(out)

}

K1_adj_T4<-function(t, mu, g, q){


	ks<-Get_Ks(mu, g)
	k1 = ks[,1]
	k2 = ks[,2]
	k3 = ks[,3]
	k4 = ks[,4]

	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp<-k1 + 2* k2 * t1/2 + 3* k3* t1^2/6 + 4* k4 *t1^3 / 24 
		out[i]<-sum(temp) -q
	}
	return(out)

}

K2_T4<-function(t, mu, g){


	ks<-Get_Ks(mu, g)
	k1 = ks[,1]
	k2 = ks[,2]
	k3 = ks[,3]
	k4 = ks[,4]

	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp<-2* k2 /2 + 3*2* k3* t1/6 + 4*3* k4 *t1^2 / 24
		out[i]<-sum(temp)
		
	}
	return(out)

}




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
		temp1<-1 - mu + mu* exp(g * t1)
		temp2<-mu *g*exp(g*t1)
		out[i]<-sum(temp2/temp1) -q 
	}
	
	
	return(out)
}

K1_adj_Fast<-function(t, mu, g, q){

	A1<-exp(g %x% t(t))
	temp1<-1-mu + mu * A1
	temp2<-mu * g * A1
	out<-colSums(temp2/temp1)-q
	return(out)
}

#t<-seq(-1000,1000, length.out=10)
#out1<-K1_adj(t, mu, g, qinv)

K2<-function(t, mu, g){

	n.t<-length(t)
	out<-rep(0,n.t)
	
	for(i in 1:n.t){
		t1<-t[i]
		temp1<-(1 - mu + mu* exp(g * t1))^2
		temp2<-(1-mu) * mu * g^2 * exp(g*t1)
		out[i]<-sum(temp2/temp1)
	}
	return(out)
}

Get_Spline_Fun<-function(mu, g, n){
	
	node<-c(-2*n, -n, -sqrt(n), 0, sqrt(n), n, 2*n)
	y1<-Korg(node ,mu ,g)
	sfun<-splinefun(node , y1 ,method = "natural")

	return(sfun)
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
	
	temp1<-zeta * q - k1

	
	w<-sign(zeta) * (2 *temp1)^{1/2}
	v<- zeta * (k2)^{1/2}
			
	Z.test<-w + 1/w * log(v/w)
	
	if(Z.test > 0){
		pval<-pnorm(Z.test, lower.tail = FALSE)
	} else {
		pval<-pnorm(Z.test, lower.tail = TRUE)
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
			
	#
	if(abs(var1) < 10^-5){	
		
		pval=0.5
	
	} else if(abs(q - m1)/sqrt(var1) < Cutoff || abs(-q - m1)/sqrt(var1) < Cutoff){
		
		pval=pval.noadj
		
	} else {
		
		out.uni1<-try(uniroot(K1_adj, c(-100,100), mu=mu, g=g, q=q), silent=TRUE)
		out.uni2<-try(uniroot(K1_adj, c(-100,100), mu=mu, g=g, q=qinv), silent=TRUE)
		
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
			
			if(pval < 10^-3 && pval.noadj > 0.1){
				pval<-pval.noadj
				Is.converge=FALSE				
			}
		
		} else {
			pval<-pval.noadj
			Is.converge=FALSE
		}
		
		
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
	idx<-which(g> 0)
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
		
		out<-uniroot(K1_adj_S, c(-100,100), q=q, sfun=sfun)
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
	} 

	
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

ScoreTest_wSaddleApprox<-function(G, obj.null, Cutoff=2){

	#Cutoff=2
	if(Cutoff < 10^-2){
		Cutoff=10^-2
	}
	
	G1<-G - mean(G)
	n.g<-sum(G)
	Y<-obj.null$y
	pi<-obj.null$fitted.values

	q<-sum(G1 * (Y-pi)) /sqrt(n.g)
	out<-Saddle_Prob(q, mu=pi, g=G1/sqrt(n.g), Cutoff=Cutoff)
	return(out)
}


##################################################################
## MAIN FUNCTION:  single.b.firth
##################################################################
single.b.spa <- function() {
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
        obj.null <- glm(pheno ~ cov, family=binomial)
        
        for (i in 1:m) {
            re <- ScoreTest_wSaddleApprox(genos[i,,drop=FALSE], obj.null, Cutoff=2)
            p[i] <- re$p.value
            
            add[i,1:2] <- c(re$p.value.NA, re$Is.converge)
        }
    }

    return(list(p=p, add=add, cname=cname))
}


