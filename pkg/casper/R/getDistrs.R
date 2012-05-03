
firstBamReads <- function(bam, nreads) {
  if (nreads < length(bam$qname)) {
    id <- unique(bam$qname[1:nreads])
    sel <- bam$qname %in% id
    bam <- lapply(bam,function(z) z[sel])
  }
  return(bam)
}

getDistrs<-function(DB, bam, nreads=4*10^6){
  
  #Format exons as RangedData
  cat("Calculating fragment length distribution\n")

  exonsRD <- DB@exonsNI
  
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
  sel <- unlist(lapply(DB@transcripts, length))==1
  oneTx <- DB@transcripts[sel]
  oneTx <- unlist(oneTx, recursive=F)
  names(oneTx) <- sub("[0-9]+\\.", "", names(oneTx))
  
  oneExons <- exonsRD[exonsRD$id %in% unlist(oneTx),]
  oneExons <- as.data.frame(oneExons)
  rownames(oneExons) <- oneExons$id
  oneExons <- oneExons[as.character(unlist(oneTx)),]
  islandStrand <- DB@islandStrand

  exon_rank <- unlist(lapply(oneTx, function(x) 1:length(x)))
  exon_strand <- islandStrand[as.character(DB@exon2island$island[match(oneExons$id, DB@exon2island$id)])]
  oneExons$strand <- exon_strand
  oneExons$exon_rank <- exon_rank

  txid <- cumsum( oneExons$exon_rank==1 )
  signwidth <- ifelse(oneExons$strand=="+", 1, -1)
  signwidth <- signwidth*(oneExons$end - oneExons$start + 1)
  oneExons$exstnogap <- 1+unlist(tapply(signwidth,INDEX=txid,function(z) {if(!any(z<0)){ c(0,cumsum(z[-length(z)]))} else {t <- -z ; rev(c(0, cumsum(rev(t)[-length(t)])))  }}))
  oneExons$txlength <- unlist(tapply((oneExons$end - oneExons$start + 1),INDEX=txid,function(z) rep(sum(z),length(z))))

  oneExons <- RangedData(oneExons)
  over <- findOverlaps(frags, oneExons, type="within")
  
  #tx <- Exons[names(oneTx)]
  #ExonsRD <- RangedData(tx@unlistData)  
  #over1 <- findOverlaps(frags, ExonsRD, type="within")
  #txid <- cumsum(ExonsRD$exon_rank==1)
  #ExonsRD$exstnogap <- 1+unlist(tapply(width(ExonsRD),INDEX=txid,function(z) c(0,cumsum(z[-length(z)]))))
  #ExonsRD$txlength <- unlist(tapply(width(ExonsRD),INDEX=txid,function(z) rep(sum(z),length(z))))
  
  ##to be removed
  #exstnogap1 <- ExonsRD$exstnogap[subjectHits(over1)]
  #txlength1 <- ExonsRD$txlength[subjectHits(over1)]
  #exst1 <- start(ExonsRD)[subjectHits(over1)]
  #exen1 <- end(ExonsRD)[subjectHits(over1)]
  #readst1 <- start(frags)[queryHits(over1)]
  #readen1 <- end(frags)[queryHits(over1)]
  #str1 <- as.character(ExonsRD$strand[subjectHits(over1)])
  
  exstnogap <- oneExons$exstnogap[subjectHits(over)]
  txlength <- oneExons$txlength[subjectHits(over)]
  exst <- start(oneExons)[subjectHits(over)]
  exen <- end(oneExons)[subjectHits(over)]
  readst <- start(frags)[queryHits(over)]
  readen <- end(frags)[queryHits(over)]
  str <- oneExons$strand[subjectHits(over)]

  stDis <- double(length(readst))
  sel <- str=='+';  stDis[sel] <- (exstnogap[sel]+readst[sel]-exst[sel])/txlength[sel]
  sel <- str=='-'; stDis[sel] <- (exstnogap[sel]+exen[sel]-readen[sel])/txlength[sel]

  #sel <- str=='+';
  #stDisP <- double(sum(sel)); stDisP <- (exstnogap[sel]+readst[sel]-exst[sel])/txlength[sel]
  #sel <- str=='-';
  #stDisM <- double(sum(sel)); stDisM <- (exstnogap[sel]+exen[sel]-readen[sel])/txlength[sel]
  #sel <- str1=='+';
  #stDisP1 <- double(sum(sel)); stDisP1 <- (exstnogap1[sel]+readst1[sel]-exst1[sel])/txlength1[sel]
  #sel <- str1=='-';
  #stDisM1 <- double(sum(sel)); stDisM1 <- (exstnogap1[sel]+exen1[sel]-readen1[sel])/txlength1[sel]
  
  
  #x11()
  #plot(density(stDis), ylim=c(0,2))
  #lines(density(distrs$stDis), col="red")
  #lines(density(stDisM), col="red")
  #plot(density(stDisP1), ylim=c(0,2), col="green")
  #lines(density(stDisM1), col="blue")

  
  list(lenDis=ld, stDis=stDis)
}

