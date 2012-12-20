require(methods)

buildRD<-function(reads){
  if(sum(grepl("chr", unique(reads$chrom)))>0) {
    reads<-GRanges(ranges=IRanges(start=ifelse(reads$end<reads$start, reads$end, reads$start), end=ifelse(reads$end>reads$start, reads$end, reads$start)), seqnames=reads$chrom, id=reads$key, flag=reads$flag, rid=reads$rid, strand=reads$rstrand, XS=reads$strand)
  } else {
    reads<-GRanges(ranges=IRanges(start=ifelse(reads$end<reads$start, reads$end, reads$start), end=ifelse(reads$end>reads$start, reads$end, reads$start)), seqnames=paste("chr", reads$chrom, sep=""), id=reads$key, flag=reads$flag, rid=reads$rid, strand=reads$rstrand, XS=reads$strand)
  }
  res<-reads
  res
}

uniquifyQname<-function(bam, seed=1){   	
  qname <- vector(mode='character', length=length(bam$qname))
  dups <- .Call("uniqQname", bam$qname, length(bam$qname), bam$pos, bam$mpos, qname)
  bam$qname <- dups[[1]]
  if(length(dups[[2]])>0){
    probPos<-bam$qname %in% dups[[2]]
    probID<-paste(bam$qname[probPos], bam$pos[probPos])
    set.seed(seed)
    sel<-unlist(tapply(1:sum(probPos), probID, function(x) x[sample(1:length(x), 1)]))
    fixed<-lapply(lapply(bam, "[", probPos), "[", sel)
    bam1<-lapply(bam, function(x) x[!probPos])    
    res<-lapply(names(bam), function(x) c(bam1[[x]], fixed[[x]]))
    names(res)<-names(bam)
  } else {
    res<-bam
  }
  qdup<-bam$qname[bam$qname %in% bam$qname[duplicated(bam$qname)]] 
  uni<-bam$qname %in% qdup
  res<-lapply(res, "[", uni)
  names(res)<-names(bam)
  res
}

nbReads <- function(bam0) {
    tab <- table(bam0$cigar)
      count <- sapply(gregexpr("M",names(tab)),length)
      tjunx <- sum(tab*(count-1))
      treads <- sum(tab*count)
      c(tjunx=tjunx, treads=treads)
  }

procB <- function(bam, strnd, seed=1){
  
  lev <- levels(bam$strand)
  bam$rname<-as.character(bam$rname)
  bam<-uniquifyQname(bam, seed)
  cat("Calculating total number of reads...\n")
  
  nreads <- nbReads(bam)
  njunx <- nreads['tjunx']
  nreads <- nreads['treads']+10

  cat("done.\nProcessing cigars and building read's object...\n")
  len=vector(mode="integer", length=nreads)
  flag=vector(mode="integer", length=nreads)
  strs=vector(mode="integer", length=nreads)
  rid=vector(mode="integer", length=nreads)
  key=vector(mode="character", length=nreads)
  chrom=vector(mode="character", length=nreads)
  strand=vector(mode="integer", length=nreads)
  jchrom=vector(mode="character", length=njunx)
  jstrs=vector(mode="integer", length=njunx)
  jlen=vector(mode="integer", length=njunx)
  data<-.Call("procBam", bam$qname, bam$flag,  bam$rname, bam$pos, bam$cigar, as.numeric(bam$strand), length(bam$pos), nreads, njunx, len, strs, flag, key, chrom, rid, strand, jchrom, jstrs, jlen)
  names(data) <- c('end', 'start', 'flag', 'key', 'chrom', 'rid', 'rstrand', 'jchrom', 'jstart', 'jend')
  junx <- GRanges(ranges=IRanges(start=data$jstart[data$jstart!=0], end=data$jend[data$jstart!=0]), seqnames=data$jchrom[data$jstart!=0], XS=rep(strnd, sum(data$jstart!=0)))
  data[[7]] <- lev[data[[7]]]
  sel<-data[[1]]!=0 & !is.na(data[[3]])
  data<-lapply(data, "[", sel)
  data[[length(data)+1]] <- rep(strnd, length(data[[1]]))
  names(data)[length(data)] <- 'strand'
  cat("done.\n")
  ans<-buildRD(data)
  id <- values(ans)$id
  idnum <- rep(0,length(id))
  idnum[which(id[-length(id)] != id[-1])+1] <- 1
  idnum <- cumsum(idnum)
  values(ans)$id <- idnum
  list(pbam=ans, junx=junx)
}

procBam<-function(bam, stranded=FALSE, seed=1){
  require(IRanges)
  if(stranded) {
    minus <- bam$tag$XS=='-'
    mjunx <- GRanges(IRanges(0,0), space=NULL)
    if(sum(minus)>0){
      minus <- lapply(bam, '[', minus)
      minus <- procB(minus, "-", seed=seed)
      mjunx <- minus$junx
      minus <- minus$pbam
    } else minus <- GRanges(IRanges(0,0), space=NULL)
    plus <- bam$tag$XS=='+'
    mplus <- GRanges(IRanges(0,0), space=NULL)
    if(sum(plus)>0){
      plus <- lapply(bam, '[', plus)
      plus <- procB(plus, "+", seed=seed)
      pjunx <- plus$junx
      plus <- plus$pbam
    } else plus <- GRanges(IRanges(0,0), space=NULL)
    ans <- list(plus=plus, pjunx=pjunx, minus=minus, mjunx=mjunx, stranded=TRUE)
  } else {
    pbam <- procB(bam, "*", seed=seed)
    ans <- list(pbam=pbam$pbam, junx=pbam$junx, stranded=FALSE)
  }
  if(!ans$stranded) ans <- new("procBam",pbam=ans$pbam, junx=ans$junx,stranded=ans$stranded)
  else ans <- new("procBam",pbam=list(plus=ans$plus, minus=ans$minus),junx=list(plus=ans$pjunx, minus=ans$mjunx), stranded=ans$stranded)
}
