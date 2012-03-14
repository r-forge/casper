calcDenovo <- function(distrs, genomeDB, pc, readLength, geneid, priorq=3, minpp=0.01, verbose=FALSE) {
  #to do: based function on islands rather than genes
  if (missing(readLength)) stop("readLength must be specified")
  #if (class(genomeDB)!='') stop("class of genomeDB must be...")
  startcdf <- as.double(ecdf(distrs$stDis)(seq(0,1,.001)))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  minpp <- as.double(minpp)
  verbose <- as.integer(verbose)
  priorprob <- function(nexonsGene, nexonsVariant) { return(1) }
  if (!missing(geneid)) {
    args <- getDenovoArgs(genomeDB=genomeDB,pc=pc,geneid=geneid)
    ans <- calcDenovoSingle(exons=args$exons,exonwidth=args$exonwidth,transcripts=args$transcripts,geneid=as.integer(geneid),pc=args$pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorprob=priorprob,priorq=priorq,minpp=minpp,verbose=verbose)
  } else {
  }
  colnames(ans[[1]]) <- c('model','posprob')
  ans[[2]] <- data.frame(ans[[2]],ans[[3]])
  ans[[3]] <- NULL
  colnames(ans[[2]]) <- c('model','expr','varName')
  ans[[3]] <- data.frame(exon=ans[[3]],varName=ans[[4]])
  ans[[4]] <- NULL
  names(ans) <- c('posprob','expression','variants')
  return(ans)
}


calcDenovoSingle <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, priorprob, priorq, minpp, verbose){
  ans <- .Call("calcDenovo",as.integer(exons),as.integer(exonwidth),transcripts,geneid,pc,startcdf,lendis,lenvals,readLength,priorprob,priorq,minpp,verbose)
  return(ans)
}

getDenovoArgs <- function(genomeDB, pc, geneid) {
#Returns arguments needed by calcDenovo
  #Select known transcripts for given gene
  sel <- unlist(values(genomeDB$txs)[['gene_id']]==geneid)
  if (any(sel)) {
    txName <- unlist(values(genomeDB$txs)[['tx_name']])[sel]
    transcripts <- genomeDB$newTx[txName]
    #Select list of exon ids & widths for given gene
    exons <- unique(unlist(transcripts)); exons <- exons[order(exons)]  #to do: order based on strand. take exons from genomeDB$gene2exon[[geneid]]
    sel <- genomeDB$exonsNI[['id']] %in% exons
    exonwidth <- as.integer(width(genomeDB$exonsNI)[sel])
    #Select path counts for given gene
    n <- strsplit(names(pc),split='-|\\.')
    n <- lapply(n,'[',-1)
    exonchar <- as.character(exons)
    sel <- sapply(n,function(z) any(z %in% exonchar))
    #Type coercion
    exons <- as.integer(exons)
    exonwidth <- as.integer(exonwidth)
    transcripts <- lapply(transcripts,as.integer)
    pc <- pc[sel]; n <- names(pc)
    pc <- as.integer(pc); names(pc) <- n
    ans <- list(exons=exons,exonwidth=exonwidth,transcripts=transcripts,pc=pc)
  } else {
    ans <- list(exons=NULL,exonwidth=NULL,transcripts=NULL,pc=NULL)
  }
  return(ans)
}

