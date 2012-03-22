calcDenovo <- function(distrs, genomeDB, pc, readLength, geneid, priorq=3, minpp=0.01, selectBest=FALSE, verbose=FALSE, mc.cores=1) {
  if (missing(readLength)) stop("readLength must be specified")
  if (class(genomeDB)!='denovoGenome') stop("genomeDB must be of class 'denovoGenome'")
  startcdf <- as.double(ecdf(distrs$stDis)(seq(0,1,.001)))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  minpp <- as.double(minpp)
  selectBest <- as.integer(selectBest)
  verbose <- as.integer(verbose)
  priorprob <- function(nexonsGene, nexonsVariant) { return(1) }
  if (missing(geneid)) geneid <- names(genomeDB@genes)[sapply(genomeDB@genes,length)>1]
  f <- function(z) {
    genesel <- genomeDB@genes[[z]]
    exons <- as.integer(names(genesel))
    pc <- pc[[z]]
    ans <- calcDenovoSingle(exons=exons,exonwidth=width(genesel),transcripts=genomeDB@transcripts[[z]],geneid=as.integer(z),pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorprob=priorprob,priorq=priorq,minpp=minpp,selectBest=selectBest,verbose=verbose)
    formatDenovoOut(ans,genesel)
  }

  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      ans <- mclapply(geneid,f,mc.cores=mc.cores)
    } else stop('multicore library has not been loaded!')
  } else {
    ans <- lapply(geneid,f)
  }
  return(ans)
}

formatDenovoOut <- function(ans, genesel) {
  colnames(ans[[1]]) <- c('model','posprob')
  ans[[2]] <- data.frame(ans[[2]],ans[[3]])
  ans[[3]] <- NULL
  colnames(ans[[2]]) <- c('model','expr','varName')
  variants <- genesel[as.character(ans[[3]])]
  names(variants) <- NULL
  ans[[3]] <- RangedData(ranges=variants,space=ans[[4]])
  ans[[4]] <- NULL
  names(ans) <- c('posprob','expression','variants')
  return(ans)
}

calcDenovoSingle <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, priorprob, priorq, minpp, selectBest, verbose) {
  ans <- .Call("calcDenovo",as.integer(exons),as.integer(exonwidth),transcripts,geneid,pc,startcdf,lendis,lenvals,readLength,priorprob,priorq,minpp,selectBest,verbose)
  return(ans)
}



