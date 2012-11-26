require(methods)

buildRD<-function(reads){
    if(sum(grepl("chr", unique(reads[[5]])))>0) {
        reads<-RangedData(IRanges(start=ifelse(reads[[1]]<reads[[2]], reads[[1]], reads[[2]]), end=ifelse(reads[[1]]>reads[[2]], reads[[1]], reads[[2]])), space=reads[[5]], id=reads[[4]], flag=reads[[3]], rid=reads[[6]], strand=reads[[7]], XS=reads[[8]])
    } else {
        reads<-RangedData(IRanges(start=ifelse(reads[[1]]<reads[[2]], reads[[1]], reads[[2]]), end=ifelse(reads[[1]]>reads[[2]], reads[[1]], reads[[2]])), space=paste("chr", reads[[5]], sep=""), id=reads[[4]], flag=reads[[3]], rid=reads[[6]], strand=reads[[7]], XS=reads[[8]])
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
    sum(tab*count)
}

procB <- function(bam, strnd, seed=1){
  lev <- levels(bam$strand)
  bam$rname<-as.character(bam$rname)
  bam<-uniquifyQname(bam, seed)
  cat("Calculating total number of reads...\n")
  nreads<-nbReads(bam)+10
  cat("done.\nProcessing cigars and building read's object...\n")
  len=vector(mode="integer", length=nreads)
  flag=vector(mode="integer", length=nreads)
  strs=vector(mode="integer", length=nreads)
  rid=vector(mode="integer", length=nreads)
  key=vector(mode="character", length=nreads)
  chrom=vector(mode="character", length=nreads)
  strand=vector(mode="integer", length=nreads)
  data<-.Call("procBam", bam$qname, bam$flag,  bam$rname, bam$pos, bam$cigar, as.numeric(bam$strand), length(bam$pos), nreads, len, strs, flag, key, chrom, rid, strand)
  data[[7]] <- lev[data[[7]]]
  sel<-data[[1]]!=0 & !is.na(data[[3]])
  data<-lapply(data, "[", sel)
  data[[length(data)+1]] <- rep(strnd, length(data[[1]]))
  cat("done.\n")
  ans<-buildRD(data)
  id <- ans$id
  idnum <- rep(0,length(id))
  idnum[which(id[-length(id)] != id[-1])+1] <- 1
  idnum <- cumsum(idnum)
  ans$id <- idnum
  ans
}

procBam<-function(bam, stranded=FALSE, seed=1){
  require(IRanges)
  if(stranded) {
    minus <- bam$tag$XS=='-'
    if(sum(minus)>0){
      minus <- lapply(bam, '[', minus)
      minus <- procB(minus, "-", seed=seed)
    } else minus <- RangedData(IRanges(0,0), space=NULL)
    plus <- bam$tag$XS=='+'
    if(sum(plus)>0){
      plus <- lapply(bam, '[', plus)
      plus <- procB(plus, "+", seed=seed)
    } else plus <- RangedData(IRanges(0,0), space=NULL)
    ans <- list(plus=plus, minus=minus, stranded=TRUE)
  } else {
    pbam <- procB(bam, "*", seed=seed)
    ans <- list(pbam=pbam, stranded=FALSE)
  }
  new("procBam",pbam=ans$pbam,stranded=ans$stranded)
}
