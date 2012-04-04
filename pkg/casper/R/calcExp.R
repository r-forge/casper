
calcExp <- function(distrs, genomeDB, pc, readLength, geneid, priorq=3, mc.cores=1) {
  if (missing(readLength)) stop("readLength must be specified")
#  if (class(genomeDB)!='knownGenome') stop("genomeDB must be of class 'knownGenome'")
  
  #Format input
  startcdf <- as.double(ecdf(distrs$stDis)(seq(0,1,.001)))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  if (missing(geneid)) geneid <- names(genomeDB@genes)[sapply(genomeDB@genes,length)>1]

  exons <- lapply(genomeDB@genes,function(z) as.integer(names(z)))
  exonwidth <- lapply(genomeDB@genes,width)

  if (!all(geneid %in% names(exons))) stop('geneid not found in genomeDB@exons')
  if (!all(geneid %in% names(pc))) stop('geneid not found in pc')
  if (!all(geneid %in% names(genomeDB@transcripts))) stop('geneid not found in genomeDB@transcripts')

  #Define basic function
  f <- function(z) {
    geneid <- as.integer(z)
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- genomeDB@transcripts[z]
    pc <- pc[z]
    ans <- calcKnownMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,geneid=geneid,pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorq=priorq)
    ans
  }

  #Run
  if (mc.cores>1) {
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
  ans
}


calcKnownMultiple <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, priorq) {
  ans <- .Call("calcKnownMultiple",exons,exonwidth,transcripts,geneid,pc,startcdf, lendis, lenvals, readLength, priorq)
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
