### Wrapper R script to perform single variant analysis
is.binary <- function(x) {
  ux <- unique(x)
  for(u in ux) {
    if ( ( !is.na(u) ) && ( u != 0 ) && ( u != 1 ) ) {
      return(FALSE)
    }
  }
  return(TRUE)
}

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
if ( test == "single.q.emmaxR" ) {  ## special treatment for emmaxR
  eigf <- args[18]
  remlf <- args[19]
} else if ( test == "single.b.spagmmat" ) {
  varRatio <- args[18]
  obj.SA.null.file <- args[19]
  formula <- args[20]
}

setwd(paste(bindir,"/share/EPACTS",sep=""))

### Wrapper R script to perform single variant analysis
if ( !require(epactsR,lib.loc=paste(bindir,"/lib/",sep="") ) ) {
  stop(paste("Cannot find epactsR package at ",bindir,"/lib/",sep=""))
}
## the DLL file has to be loaded explictly for some machine - don't know why
dyn.load(paste(bindir,"/lib/epactsR/libs/epactsR",.Platform$dynlib.ext, sep=""))
  
pheno <- as.double(read.table(phenof)[,2])
pheno <- pheno - min(pheno,na.rm=T)
if ( covf == "NULL" ) {
  cov <- NULL;
} else {
  cov <- as.matrix(read.table(covf)[,-1])
}

binaryFlag <- is.binary(pheno)

ind <- as.integer(read.table(indf)[,2])-1
#indids <- as.character(read.table(indf)[,1])

G <- .Call("readVcf",vcf,region,field,passOnly,ind,NULL)
rnames <- rownames(G)
rownames(G) <- seq_along(rnames)
m <- nrow(G)
if ( m > 0 ) {
  n <- ncol(G)
  cov <- cbind(rep(1,n),cov)  
  NS <- rowSums(!is.na(G))
  AC <- rowSums(G,na.rm=T)
  MAC <- apply(cbind(AC,NS+NS-AC),1,min)
  MAF <- MAC/NS/2
  #MAF <- 0.5 - abs(AC/NS/2 - 0.5)
  CR <- NS/n
  sqAC <- rowSums(G*G,na.rm=T)
  varAC <- sqAC/NS - AC*AC/NS/NS
  rsq <- (sqAC-AC*AC/(NS+1e-30))/(AC - AC*AC/(2*NS+1e-30))
  rsq[rsq>1] <- 1
  if ( minRSQ > 0 ) {
    vids <- which((varAC > 0) & (MAF >= minMAF) & (MAF <= maxMAF ) & ( MAC >= minMAC ) & ( MAC <= maxMAC) & (CR >= minCallRate) &  (rsq >= minRSQ))
  } else {
    vids <- which((varAC > 0) & (MAF >= minMAF) & (MAF <= maxMAF) & ( MAC >= minMAC ) & ( MAC <= maxMAC) & (CR >= minCallRate))
  }
  genos <- G[vids,,drop=FALSE]

  #print(sum(MAF >= minAF))
  #print(MAF)
  #print(minAF)
  #print(sum(NS+NS-AC>2))
  #print(sum(CR >= minCallRate))

  source(paste(test,".R",sep="",collapse=NULL))
  #print(test)

  func <- match.fun(test)
  r <- func()
  

  if ( binaryFlag ) {
    binary.cname <- c("NS.CASE","NS.CTRL","AF.CASE","AF.CTRL") #,"CNT.CASE","CNT.CTRL");
    binary.ncols <- length(binary.cname)
    binary.add <- matrix(NA,nrow(G),4)
    is.case <- (pheno == 1)
    binary.add[,1] <- rowSums(!is.na(G[,is.case,drop=FALSE]))
    binary.add[,2] <- rowSums(!is.na(G[,!is.case,drop=FALSE]))
    binary.add[,3] <- rowSums(G[,is.case,drop=FALSE],na.rm=TRUE)/binary.add[,1]/2
    binary.add[,4] <- rowSums(G[,!is.case,drop=FALSE],na.rm=TRUE)/binary.add[,2]/2
    binary.add[binary.add[,1] == 0,3] <- NA
    binary.add[binary.add[,2] == 0,4] <- NA
    #binary.add[,5] <- apply(t(apply(
    #ina <- is.na(G)
    #binary.add[,1] <- rowSums((matrix(pheno,nrow(G),ncol(G),byrow=T) == 1) * (1-ina), na.rm=T)
    #binary.add[,2] <- rowSums((matrix(pheno,nrow(G),ncol(G),byrow=T) == 0) * (1-ina), na.rm=T)
    #binary.add[,3] <- rowSums((matrix(pheno,nrow(G),ncol(G),byrow=T) == 1) * G * (1-ina), na.rm=T)/binary.add[,1]
    #binary.add[,4] <- rowSums((matrix(pheno,nrow(G),ncol(G),byrow=T) == 0) * G * (1-ina), na.rm=T)/binary.add[,2]
    #binary.add[binary.add[,1] == 0,3] <- NA
    #binary.add[binary.add[,2] == 0,4] <- NA
  } else {
    binary.cname <- NULL
    binary.ncols <- 0
  }

  out <- matrix(NA,m,5+(minRSQ>0)+ncol(r$add)+binary.ncols)
  
  out[,1] <- NS
  out[,2] <- AC
  out[,3] <- CR
  out[,4] <- MAF
  if ( minRSQ > 0 ) {
    colnames(out) <- c("NS","AC","CALLRATE","MAF","PVALUE",r$cname,"RSQ",binary.cname);
    out[,6+ncol(r$add)] <- rsq
  } else {
    colnames(out) <- c("NS","AC","CALLRATE","MAF","PVALUE",r$cname,binary.cname);
  }
  out[vids,5] <- r$p
  out[vids,5+1:ncol(r$add)] <- r$add
  if ( binaryFlag ) {
    out[,5+ncol(r$add)+1:4] <- binary.add
  }
  rownames(out) <- rnames
  
  if (class(try(rnames[seq_along(rnames)], silent = TRUE)) == "try-error") {
    tmp <- rownames(out)
    goodnames <- sapply(seq_along(tmp), function (iname) {
      class(try(tmp[iname], silent = TRUE)) != "try-error"
    })
    out <- out[which(goodnames), ]
    cat(paste(sum(!goodnames), 'variants have been excluded due to error: "Value of SET_STRING_ELT() must be a \'CHARSXP\' not a \'integer\'"'), file = paste0(outf, ".log"))
  } ##else {}
  
  print(warnings())
  .Call("writeMatrix",out,outf)
} else {
  write.table(NULL,outf,row.names=F,col.names=F)  
}
