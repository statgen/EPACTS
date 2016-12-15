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

## group.q.reverse.wilcox() : Reverse regression with wilcox test
## KEY FEATURES : 0/1 collapsing variable ~ rare variants
## TRAITS  : QUANTITATIVE (GAUSSIAN)
## RETURNS : PVALUE, BETA, SEBETA, ZSTAT
## MISSING VALUE : IMPUTED AS MAJOR ALLELES
group.q.wilcox <- function() {
  cname <- c("STAT")
  res <- lm(pheno ~ cov - 1)$residual
  m <- nrow(genos)
  if ( m > 0 ) {
    #print(genos)
    g <- as.double(colSums(genos,na.rm=T) > 0)
    if ( var(g,na.rm=T) > 0 ) {
      #print(length(pheno))
      #print(length(g))
      r <- wilcox.test(res~g)
      return(list(p=r$p.value,add=matrix(r$statistic,1,1),cname=cname))      
    }
  }
  return(list(p=NA,add=matrix(NA,1,1),cname=cname))
}
