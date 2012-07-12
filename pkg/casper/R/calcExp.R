
calcExp <- function(distrs, genomeDB, pc, readLength, geneid, relativeExpr=TRUE, priorq=3, report=0, niter=10^3, burnin=100, mc.cores=1) {
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if (pc@denovo) stop("pc must be a pathCounts object from known genome")
  
  #Format input
  startcdf <- distrs$stDis(seq(0,1,.001))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  report <- as.integer(report)
  niter <- as.integer(niter)
  burnin <- as.integer(burnin)
  if (!report %in% 0:2) stop("Argument report must be equal to 0, 1 or 2")
  if ((niter<=burnin) & report>=2) stop("Too many burnin iterations specified")
  if (missing(geneid)) geneid <- names(genomeDB@islands)[sapply(genomeDB@islands,length)>1]

  exons <- lapply(genomeDB@islands,function(z) as.integer(names(z)))
  exonwidth <- lapply(genomeDB@islands,width)
  strand <- genomeDB@islandStrand
  
  if (!all(geneid %in% names(exons))) stop('geneid not found in genomeDB@islands')
  if (!all(geneid %in% names(pc@counts))) stop('geneid not found in pc')
  if (!all(geneid %in% names(genomeDB@transcripts))) stop('geneid not found in genomeDB@transcripts')
  if (!all(geneid %in% names(genomeDB@islandStrand))) stop('geneid not found in genomeDB@islandStrand')
    
  #Define basic function
  f <- function(z) {
    geneid <- as.integer(z)
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- genomeDB@transcripts[z]
    strand <- as.list(as.integer(ifelse(strand[z]=='+', 1, -1)))
    pc <- pc@counts[z]
    ans <- calcKnownMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,geneid=as.list(geneid),pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorq=priorq, strand=strand, report=report, niter=niter, burnin=burnin)
    if (report==0) {
      ans <- lapply(ans, function(z) { res=vector("list",1); res[[1]]= z[[1]]; names(res[[1]])= z[[2]]; res })
    } else if (report==1) {
      ans <- lapply(ans, function(z) { res=vector("list",2); res[[1]]= z[[1]]; res[[2]]= z[[3]]; names(res[[1]])= names(res[[2]])= z[[2]]; res })
    } else if (report==2) {
      ans <- lapply(ans, function(z) { res=vector("list",3); res[[1]]= z[[1]]; res[[2]]= z[[3]]; res[[3]]= matrix(z[[4]],nrow=niter-burnin); names(res[[1]])= names(res[[2]])= z[[2]]; res })
    }
    ans
  }

  #Run
  if (mc.cores>1 && length(geneid)>mc.cores) {
    if ('multicore' %in% loadedNamespaces()) {
      #split into smaller jobs
      nsplit <- floor(length(geneid)/mc.cores)
      geneid <- lapply(1:mc.cores, function(z) { geneid[((z-1)*nsplit+1):min((z*nsplit),length(geneid))] })
      ans <- mclapply(geneid,f,mc.cores=mc.cores)
      ans <- do.call(c,ans)
     } else stop('multicore library has not been loaded!')
  } else {
    ans <- f(geneid)
    names(ans) <- geneid
  }
  #Format as ExpressionSet
  transcript <- unlist(lapply(ans, function(z) names(z[[1]])))
  gene <- rep(names(ans),sapply(ans,function(z) length(z[[1]])))
  fdata <- data.frame(transcript=transcript, gene=gene)
  exprsx <- matrix(unlist(lapply(ans,'[[',1)),ncol=1)
  if (report>=1) fdata$SE <- sqrt(unlist(lapply(ans,'[[',2)))
  if (report>=2) {
    q <- lapply(ans,function(z) apply(z[[3]],2,quantile,probs=c(.025,.975)))
    q <- t(do.call(cbind,q))
    fdata$ci95.low <- q[,1]; fdata$ci95.high <- q[,2]
  }
  rownames(exprsx) <- rownames(fdata) <- fdata$transcript
  fdata <- new("AnnotatedDataFrame",fdata)
  ans <- new("ExpressionSet",exprs=exprsx,featureData=fdata)

  #Return absolute expression levels
  if (!relativeExpr) {
    nreads <- rep(sapply(pc@counts[geneid],sum),sapply(genomeDB@transcripts[geneid],length))
    exprs(ans) <- exprs(ans)*nreads
    exprs(ans) <- t(t(exprs(ans))/colSums(exprs(ans)))
  }
  ans
}


calcKnownMultiple <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, priorq, strand, report, niter, burnin) {
  ans <- .Call("calcKnownMultiple",exons,exonwidth,transcripts,geneid,pc,startcdf, lendis, lenvals, readLength, priorq, strand, report, niter, burnin)
  return(ans)
}


calcExpOld<-function(distrs, genomeDB, pc, readLength){
  if (missing(readLength)) stop("readLength must be specified")
  if (class(genomeDB)!="knownGenome") stop("genomeDB must be of class 'knownGenome'")
  dexo<-as.data.frame(genomeDB@exonsNI)
  mexo<-as.matrix(dexo[,c(5,4)])
  startcdf <- as.double(ecdf(distrs$stDis)(seq(0,1,.001)))
  #stafun = ecdf(distrs$stDis)
  #fill<-c(rep(0,(min(as.numeric(names(distrs$lenDis))))), distrs$lenDis)
  #names(fill)[1:min(as.numeric(names(distrs$lenDis)))]<-0:(min(as.numeric(names(distrs$lenDis)))-1)
  lendis<- as.double(distrs$lenDis/sum(distrs$lenDis))
  exp<-.Call("calc", mexo, genomeDB@newTx, pc, startcdf, lendis, as.integer(names(lendis)), as.integer(readLength))
  names(exp)<-names(genomeDB@newTx)
  exp<-as.matrix(exp)
  featureData<-as.numeric(unlist(genomeDB@txs@elementMetadata$gene_id))
  names(featureData)<-unlist(genomeDB@txs@elementMetadata$tx_name)
  featureData<-as.data.frame(as.matrix(featureData[rownames(exp)]))
  metadata <- data.frame(labelDescription = "gene_id", row.names="gene_id")
  featureData<-new("AnnotatedDataFrame", data=featureData, varMetadata=metadata)
  exp<-new("ExpressionSet", exprs=exp, featureData=featureData, annotation="refseq")
  names(fData(exp)) <- 'entrezid'
  #Compute rpkm
  exons <- unlist(genomeDB@exons)
  geneid <- cumsum(values(exons)$exon_rank==1)
  txlength <- tapply(width(exons),geneid,FUN=sum)
  names(txlength) <- names(genomeDB@exons)
  fData(exp)$txlength <- txlength[featureNames(exp)]
  rpkm <- matrix(10^9 * exprs(exp)[,] / fData(exp)$txlength, nrow=nrow(exp), ncol=ncol(exp))
  exp<-new("ExpressionSet", exprs=rpkm, featureData=featureData(exp), annotation="refseq")
  exp
}
