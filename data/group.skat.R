##################################################################
## GENE-LEVEL BURDEN TEST
## INPUT VARIABLES: 
##   n        : total # of individuals
##   genos    : genotype matrix for each gene
##   NS       : number of called samples for each marker
##   AC       : allele count for each marker
##   MAC      : minor allele count for each marker
##   MAF      : minor allele frequency
##   vids     : indices from 1:n after AF/AC threshold
## EXPECTED OUTPUT : list(p, addcols, addnames) for each genos row
##   p        : p-value
##   add      : additional column values
##   cname    : additional column names
##################################################################

## check the installation of SKAT package
if ( !require(mmSKAT,lib.loc=paste(bindir,"/lib/",sep="") ) ) {
  stop("Cannot find mmSKAT package");
#  if ( batch == 0 ) { # install for the first time
#    install.packages('SKAT_0.77.tar.gz');
#    if(!require(SKAT)) { stop("Package SKAT is not found") }
#  }
#  else {
#    found <- FALSE;
#    for(i in 1:10) {
#      print("Waiting for SKAT package being installed..")
#      Sys.sleep(5);
#      if(require(SKAT)) { found <- TRUE; break; }
#    }
#    if ( found ) { print("Waiting for SKAT package being installed..") }
#    else { stop("Package SKAT is not found after a long wait"); }
#  }
}

library(mmSKAT)

## group.q.skat() : SKAT implementation
## KEY FEATURES : 0/1 collapsing variable ~ rare variants
## TRAITS  : QUANTITATIVE (GAUSSIAN)
## RETURNS : SKAT_PVALUE, SKAT_QSTAT, rSKAT_PVALUE, rSKAT_QSTAT, mdSKAT_PVALUE, mdSKAT_QSTAT, SKATO_PVALUE, SKATO_RHO, rSKATO_PVALUE, rSKATO_RHO for SKAT, rSKAT, SKAT with Madsen Browning type weights, SKAT_O, rSKAT_O
## MISSING VALUE : IMPUTED AS MAJOR ALLELES

## SKAT
## binary vs quantitiative
## weighted vs flat
## optimal vs original
group.skat <- function() {
  #library(SKAT)
  #print(packageDescription("SKAT"))
  if ( skatOptimal ) {
    cname <- c("STATRHO");
  }
  else {
    cname <- c("QSTAT");
  }
  
  #cname <- c("QTSTAT","rSKAT_PVALUE","rSKAT_QSTAT","mdSKAT_PVALUE","mdSKAT_QSTAT","SKATO_PVALUE","SKATO_RHO","rSKATO_PVALUE","rSKATO_RHO")
  ## flip ref/alt to minor/major
  flip <- which(AC>NS)
  genos[flip,] <- 2-genos[flip,]
  m <- nrow(genos)

  genos.nona <- genos
  genos.nona[is.na(genos)] <- 0
  if ( qr(genos.nona)$rank < 1 ) {  ## Change suggested by Jason
    m <- 0
  }
  if ( m > 0 ) {
    if ( binaryFlag ) {
      obj <- SKAT_Null_Model(pheno~cov-1,out_type="D",n.Resampling=0, type.Resampling="bootstrap", Adjustment=skatAdjust)
      # Dichotomous outcome
      #obj <- SKAT_Null_Model(pheno~cov-1,out_type="C");
    }
    else {
      obj <- SKAT_Null_Model(pheno~cov-1,out_type="C",n.Resampling=0, type.Resampling="bootstrap", Adjustment=skatAdjust) # Continuous outcome
    }
    # Default SKAT with linear weighted kernel and weights = beta(MAF,1,25)
    if ( skatOptimal ) {
      if ( skatFlat ) {
        if ( skatAdjust ) {
          r <- SKAT(t(genos),obj,kernel="linear",method="optimal.adj",impute.method="fixed")
        }
        else {
          r <- SKAT(t(genos),obj,kernel="linear",method="optimal",impute.method="fixed")
        }
      }
      else {
        if ( skatAdjust ) {
          r <- SKAT(t(genos),obj,kernel="linear.weighted",method="optimal.adj",weights.beta=betas,impute.method="fixed")
        }
        else {
          r <- SKAT(t(genos),obj,kernel="linear.weighted",method="optimal",weights.beta=betas,impute.method="fixed")
        }
      }

      p <- r$p.value
      if ( m == 1 ) {
        rho <- NA
      }
      else {
        rho <- r$param$rho_est
        if ( length(rho) > 1 ) {
          rho <- NA  ## temporary solution
        }
        if ( is.null(rho) ) {
          rho <- NA
        }        
      }
      return(list(p=p,add=c(rho),cname=cname))
    }
    else {
      if ( skatFlat ) {
        r <- SKAT(t(genos),obj,kernel="linear",method="davies",impute.method="fixed")
      }
      else {
        r <- SKAT(t(genos),obj,kernel="linear.weighted",method="davies",weights.beta=betas,impute.method="fixed")
        #stop("foo")                
      }
      Q <- r$Q
      p <- ifelse(r$p.value==0,r$param$liu_pval,r$p.value)
      rm(r)
      gc()
      return(list(p=p,add=c(Q),cname=cname))
    }
  }
  else {
    return(list(p=NA,add=rep(NA,length(cname)),cname=cname))
  }
}
