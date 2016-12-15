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

## group.q.reverse() : Reverse logistic regression
## KEY FEATURES : 0/1 collapsing variable ~ rare variants
## TRAITS  : QUANTITATIVE (GAUSSIAN)
## RETURNS : PVALUE, BETA, SEBETA, ZSTAT
## MISSING VALUE : IMPUTED AS MAJOR ALLELES
group.q.reverse <- function() {
  cname <- c("BETA","SEBETA","ZSTAT")
  m <- nrow(genos)
  if ( m > 0 ) {
    #print(genos)
    g <- as.double(colSums(genos,na.rm=T) > 0)
    sg <- sum(g)
    if ( ( sg > 0 ) && ( sg < n ) ) {
      r <- glm(g~pheno+cov-1,family=binomial)
      
      #print(summary(r))
      
      if ( ( r$converged ) && ( ! r$boundary ) ) {
        return(list(p=summary(r)$coefficients[1,4],
                    add=summary(r)$coefficients[1,1:3],
                    cname=cname))
      }
    }
  }
  return(list(p=NA,add=rep(NA,3),cname=cname))
}
