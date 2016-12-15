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

## group.b.burden() : Burden test using logistic regression
## KEY FEATURES : 0/1 collapsing variables to dichotomous traits
## TRAITS : BINARY
## RETURNS : PVALUE, BETA, SEBETA, ZSTAT
## MISSING VALUE : IMPUTED AS MAJOR ALLELES
group.b.burden <- function() {
  cname <- c("BETA","SEBETA","ZSTAT")
  m <- nrow(genos)
  if ( m > 0 ) {
    ## weight each allele based on allele frequency
    #weights <- 1/sqrt(NS[vids]*MAF[vids]*(1-MAF[vids]))
    g <- colSums(genos * matrix(1,nrow(genos),ncol(genos)),na.rm=T)
    if ( var(g,na.rm=T) > 0 ) {
      r <- glm(pheno~g+cov-1,family=binomial)
      
      if ( ( r$converged ) && ( ! r$boundary ) ) {
        return(list(p=summary(r)$coefficients[1,4],
                    add=summary(r)$coefficients[1,1:3],
                    cname=cname))
      }
    }
  }
  return(list(p=NA,add=rep(NA,3),cname=cname))
}
