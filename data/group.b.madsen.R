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

## group.b.madsen() : Madsen-Browning test
## KEY FEATURES : 0/1 collapsing variables to dichotomous traits
## TRAITS : BINARY
## RETURNS : PVALUE, BETA, SEBETA, ZSTAT
## MISSING VALUE : IMPUTED AS MAJOR ALLELES
group.b.madsen <- function() {
  cname <- c("STAT")
  m <- nrow(genos)
  if ( m > 0 ) {
    #print(genos)
    weights <- 1/sqrt(NS[vids]*MAF[vids]*(1-MAF[vids]))
    #print(weights)
    g <- colSums(genos * matrix(weights,nrow(genos),ncol(genos)),na.rm=T)
    if ( var(g,na.rm=T) > 0 ) {
      r <- wilcox.test(g~pheno)
      return(list(p=r$p.value,add=r$statistic,cname=cname))
    }
  }
  return(list(p=NA,add=rep(NA,1),cname=cname))
}
