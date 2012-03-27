
firstBamReads <- function(bam, nreads) {
  if (nreads < length(bam$qname)) {
    id <- unique(bam$qname[1:nreads])
    sel <- bam$qname %in% id
    bam <- lapply(bam,function(z) z[sel])
  }
  return(bam)
}

getDistrs<-function(txs, exons, bam, nreads=4*10^6){
  
  #Format exons as RangedData
  cat("Calculating fragment length distribution\n")
  exonsRD<-as.data.frame(exons@unlistData)
  exonsRD<-exonsRD[,-7]
  colnames(exonsRD)[1]<-"space"
  exonsRD<-RangedData(exonsRD)
  
  #Select a sample of reads
  bam <- firstBamReads(bam, nreads=nreads)

  #Find fragment length distribution for fragments aligning to exons larger than 1000 bases
  if (!all(c('qname','rname','qwidth','pos','mpos') %in% names(bam))) stop('bam must contain elements qname, rname, qwidth, pos, mpos')
  d <- bam$mpos - bam$pos
  sel <- d<0; n <- bam$qname[sel]; sp <- bam$rname[sel]; en <- bam$pos[sel]+bam$qwidth[sel]-1
  sel <- d>0; st <- bam$pos[sel]; names(st) <- bam$qname[sel]; st <- st[n]
  sel <- st<en; st <- st[sel]; en <- en[sel]; sp <- sp[sel]
  frags <- RangedData(IRanges(start=st,end=en),space=sp)
  n <- names(frags)[names(frags) %in% names(exonsRD)]
  fragsL<-frags[n]
  fragsL<-subsetByOverlaps(fragsL, exonsRD[width(exonsRD)>1000,], type="within")
  ld<-table(width(fragsL))
  ld <- ld[ld/sum(ld) > 0.0001]

  #Find fragment start distribution for fragments aligning to transcripts in genes with only one annotated tx

  cat("Calculating fragment start distribution\n")
  txs <- RangedData(ranges=ranges(txs),space=seqnames(txs),strand=strand(txs),values(txs))
  oneTx<- unlist(txs$gene_id) %in% names(which(table(unlist(txs$gene_id))==1))
  oneTx<-txs[oneTx,]
  oneTx<-subsetByOverlaps(oneTx, frags)
  exonsRD<-subsetByOverlaps(exonsRD, oneTx)
  over<-findOverlaps(frags, exonsRD, type="within")

  txid <- cumsum(exonsRD$exon_rank==1)
  exonsRD$exstnogap <- 1+unlist(tapply(width(exonsRD),INDEX=txid,function(z) c(0,cumsum(z[-length(z)]))))
  exonsRD$txlength <- unlist(tapply(width(exonsRD),INDEX=txid,function(z) rep(sum(z),length(z))))
   
  exstnogap <- exonsRD$exstnogap[subjectHits(over)]
  txlength <- exonsRD$txlength[subjectHits(over)]
  exst <- start(exonsRD)[subjectHits(over)]
  exen <- end(exonsRD)[subjectHits(over)]
  readst <- start(frags)[queryHits(over)]
  readen <- end(frags)[queryHits(over)]
  str <- exonsRD$strand[subjectHits(over)]
  stDis <- double(length(readst))
  sel <- str=='+'; stDis[sel] <- (exstnogap[sel]+readst[sel]-exst[sel])/txlength[sel]
  sel <- str=='-'; stDis[sel] <- (exstnogap[sel]+exen[sel]-readen[sel])/txlength[sel]
  
  list(lenDis=ld, stDis=stDis)
}
