### Wrapper R script to perform single variant analysis
#if ( !require(epactsR,lib.loc=paste(bindir,"../lib/",sep="") ) ) {
#if ( !require(epactsR) ) {
#  stop("Cannot find epactsR package");
#  found <- FALSE;
#  for(i in 1:10) {
#    install.packages('epactsR_3.1.0.tar.gz');
#    if(require(epactsR)) { found <- TRUE; break; }
#    print("Waiting for epactsR package being installed..")
#    Sys.sleep(runif(1,0,5));
#  }
#}

pheno <- as.double(read.table(phenof)[,2])
pheno <- pheno - min(pheno,na.rm=T)
n <- length(pheno)

if ( is.null(covf) || ( covf == "NULL" ) ) {
  cov <- NULL;
} else {
  cov <- as.matrix(read.table(covf)[,-1])
}
cov <- cbind(rep(1,n),cov)  
ind <- as.integer(read.table(indf)[,2])-1

source(paste(test,".R",sep="",collapse=NULL))
nr <- length(groups)

func <- match.fun(test)

#if ( test == "group.mmskat" ) {  ## special treatment for emmaxR
#  eigf <- args[18]
#  remlf <- args[19]
#}

#if ( is.null(kinf) ) {
#  K <- NULL
#} else {  ## read kinship matrix if specified (e.g. for emmaxSKAT)
#  ids <- as.character(read.table(indf)[,1])
#  K <- .Call("readKinWithIDs",kinf,ids)
#}

for(i in 1:nr) {
  ## load the genotype matrix from groupwise info
  G <- .Call("readVcfGroup",vcf,field,ind,groups[[i]],sepchr)
  NS <- rowSums(!is.na(G))
  AC <- rowSums(G,na.rm=T)
  if ( flip ) {  ## flip major/minor, AC == MAC
    fids <- which(AC>NS)
    G[fids,] <- 2-G[fids,]
    AC[fids] <- NS[fids]*2 - AC[fids]
    MAC <- AC
  } else {
    MAC <- NS - abs(NS-AC) # if AC < NS, NS-NS+AC=AC, if AC>NS, NS-AC+NS=2NS-AC
  }

  MAF <- MAC / (NS + NS)
  vids <- which((MAF >= minAF) & ( MAF <= maxAF ) & ( MAC >= minAC ))
  genos <- G[vids,,drop=FALSE]  ## informative genotypes
  if (length(vids) > 0 ) {
    collapse <- as.double(colSums((genos > 0.5),na.rm=T) > 0)
    sum.collapse <- sum(collapse,na.rm=T)
    nona.collapse <- sum(!is.na(collapse),na.rm=T)
    if ( nona.collapse > 0 ) {
      maf.collapse <- sum.collapse/nona.collapse
    } else {
      maf.collapse <- NA
    }
  } else {
    nona.collapse <- NA
    maf.collapse <- NA
  }

  #print(genos)
  #print(rowSums(genos,na.rm=T))
  r <- func()
  if ( i == 1 ) {
    out <- matrix(NA,nr,length(r$cname)+6)
    colnames(out) <- c("NS","FRAC_WITH_RARE","NUM_ALL_VARS","NUM_PASS_VARS","NUM_SING_VARS","PVALUE",r$cname)
  }
  #print(r)
  out[i,] <- c(nona.collapse,maf.collapse,nrow(G),nrow(genos),sum(AC==1),r$p,r$add)
}
rownames(out) <- paste(regions,gnames,sep="_")
print(warnings())

.Call("writeMatrix",out,outf)
