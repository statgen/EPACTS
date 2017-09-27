### Load arguments from command line
args <- commandArgs(trailingOnly=TRUE)
bindir <- args[2]  ## Directory containing the R scripts
phenof <- args[3]  ## phenotype file
covf <- args[4]    ## covariate file
indf <- args[5]    ## individual IDs (in order from VCF)
vcf <- args[6]     ## VCF file - must be bgzipped and tabixed
region <- args[7]  ## Genomic region to load
outf <- args[8]    ## Output file prefix
field <- args[9]   ## Field to extract from VCF
minMAF <- as.double(args[10])    ## Minimum MAF 
maxMAF <- as.double(args[11])    ## Maximum MAF
minMAC <- as.integer(args[12])   ## Minimum MAC
maxMAC <- as.integer(args[13])   ## Maximum MAC
minCallRate <- args[14]          ## Minimum callRate
minRSQ <- args[15]               ## 
passOnly <- (args[16] == "TRUE")
test <- args[17]

setwd(paste(bindir,"/share/EPACTS",sep=""))

### Wrapper R script to perform single variant analysis
if ( !require(epactsR,lib.loc=paste(bindir,"/lib/",sep="") ) ) {
  stop(paste("Cannot find epactsR package at ",bindir,"/lib/",sep=""))
}
## the DLL file has to be loaded explictly for some machine - don't know why
dyn.load(paste(bindir,"/lib/epactsR/libs/epactsR",.Platform$dynlib.ext, sep=""))

pnames <- scan(phenof,what=character(),nlines=1)[-1]
phenos <- as.matrix(read.table(phenof)[,-1])
cov <- as.matrix(read.table(covf)[,-1])

ind <- as.integer(read.table(indf)[,2])-1

G <- .Call("readVcf",vcf,region,field,passOnly,ind,NULL)
m <- nrow(G)
if ( m > 0 ) {
  n <- ncol(G)
  cov <- cbind(rep(1,n),cov)  
  NS <- rowSums(!is.na(G))
  AC <- rowSums(G,na.rm=T)
  MAC <- apply(cbind(AC,NS+NS-AC),1,min)
  MAF <- MAC/NS/2
  CR <- NS/n
  sqAC <- rowSums(G*G,na.rm=T)
  varAC <- sqAC/NS - AC*AC/NS/NS
  rsq <- (sqAC-AC*AC/(NS+1e-30))/(AC - AC*AC/(2*NS+1e-30))
  rsq[rsq>1] <- 1
  if ( minRSQ > 0 ) {
    vids <- which((varAC > 0) & (MAF >= minMAF) & (MAF <= maxMAF ) & ( MAC >= minMAC ) & ( MAC <= maxMAC) & (CR >= minCallRate) &  (rsq >= minRSQ))
  }
  else {
    vids <- which((varAC > 0) & (MAF >= minMAF) & (MAF <= maxMAF) & ( MAC >= minMAC ) & ( MAC <= maxMAC) & (CR >= minCallRate))
  }
  genos <- G[vids,,drop=FALSE]

  #print("foo")
  source(paste(test,".R",sep="",collapse=NULL))

  func <- match.fun(test)
  r <- func()

  print("bar")

  out <- matrix(NA,m,4+ncol(r$p)+ncol(r$add))
  ##out <- matrix(NA,m,4+2*ncol(phenos))
  out[vids,4+(1:ncol(r$p))] <- r$p
  ##out[vids,3+(1:ncol(phenos))*2] <- r$p
  out[vids,4+ncol(r$p)+(1:ncol(r$add))] <- r$add
  ##out[vids,4+(1:ncol(phenos))*2] <- r$b
  out[,1] <- NS
  out[,2] <- AC
  out[,3] <- CR
  out[,4] <- MAF
  #colnames(out) <- c("NS","AC","CR","MAF",rep(NA,2*ncol(phenos)))
  #colnames(out)[3+(1:ncol(phenos))*2] <- paste(pnames,"P",sep=".")
  #colnames(out)[4+(1:ncol(phenos))*2] <- paste(pnames,"B",sep=".")
  colnames(out) <- c("NS","AC","CR","MAF",rep(NA,ncol(r$p)),rep(NA,ncol(r$add)))
  colnames(out)[4+(1:ncol(r$p))] <- paste(pnames,"P",sep=".")
  colnames(out)[4+ncol(r$p)+(1:ncol(r$add))] <- r$cname

  rownames(out) <- rownames(G)
  print(warnings())
  .Call("writeMatrix",out,outf)
} else {
  write.table(NULL,outf,row.names=F,col.names=F)  
}
