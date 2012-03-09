calcDenovo <- function(distrs, genomeDB, pc, readLength, geneid, verbose=FALSE){
  #to do: based function on islands rather than genes
  if (missing(readLength)) stop("readLength must be specified")
  #if (class(genomeDB)!='') stop("class of genomeDB must be...")
  startcdf <- as.double(ecdf(distrs$stDis)(seq(0,1,.001)))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  verbose <- as.integer(verbose)
  priorprob <- function(nexonsGene, nexonsVariant) { return(1) }
  if (!missing(geneid)) {
    args <- getDenovoArgs(genomeDB=genomeDB,pc=pc,geneid=geneid)
    ans <- calcDenovoSingle(exons=args$exons,exonwidth=args$exonwidth,transcripts=args$transcripts,geneid=as.integer(geneid),pc=args$pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorprob=priorprob,verbose=verbose)
  } else {
  }
  return(ans)
}


calcDenovoSingle <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, priorprob, verbose){
  ans <- .Call("calcDenovo", as.integer(exons), as.integer(exonwidth), transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, priorprob, verbose)
  return(ans)
}

getDenovoArgs <- function(genomeDB, pc, geneid) {
#Returns arguments needed by calcDenovo
  #Select known transcripts for given gene
  sel <- unlist(values(genomeDB$txs)[['gene_id']]==geneid)
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
  pc <- as.integer(pc[sel])
  list(exons=exons,exonwidth=exonwidth,transcripts=transcripts,pc=pc)
}

