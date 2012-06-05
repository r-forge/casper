startDist <- function(st,fragLength,txLength) {
                                        # Estimate relative start distribution under left­truncation (st < 1 ­ fragLength/txLength)
                                        # ­ st: relative start (i.e. start/txLength)
                                        # ­ fragLength: fragment length
                                        # ­ txLength: transcript length
                                        # Output: cumulative probability function (actually, a linear interpolation)
  require(survival)
  trunc <- 1-fragLength/txLength
  sel <- trunc>st
  trunc <- trunc[sel]
  st <- st[sel]
  fit <- summary(survfit(Surv(time=1-trunc,time2=1-st,event=rep( TRUE,length(st))) ~ 1))
  s <- 1-fit$time
  pcum <- fit$surv
  f <- approxfun(s,pcum)
  sseq <- seq(0,1,.001)
  startcdf <- f(sseq)
  startcdf[1] <- 0; startcdf[length(startcdf)] <- 1
  f <- approxfun(sseq[!is.na(startcdf)], startcdf[!is.na(startcdf)])
  return(f)
}

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
  
  exstnogap <- oneExons$exstnogap[subjectHits(over)]
  txlength <- oneExons$txlength[subjectHits(over)]
  exst <- start(oneExons)[subjectHits(over)]
  exen <- end(oneExons)[subjectHits(over)]
  readst <- start(frags)[queryHits(over)]
  readen <- end(frags)[queryHits(over)]
  str <- oneExons$strand[subjectHits(over)]
  frlen <- width(frags)[queryHits(over)]
  
  stDis <- double(length(readst))
  sel <- str=='+';  stDis[sel] <- (exstnogap[sel]+readst[sel]-exst[sel])/txlength[sel]
  sel <- str=='-'; stDis[sel] <- (exstnogap[sel]+exen[sel]-readen[sel])/txlength[sel]

  stDis <- startDist(stDis, frlen, txlength)
  
  list(lenDis=ld, stDis=stDis)
}

