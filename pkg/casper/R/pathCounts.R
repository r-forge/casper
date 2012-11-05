
procPaths <- function(reads, DB, mc.cores){
    cat("Finding overlaps between reads and exons\n")
    over<-findOverlaps(reads, DB@exonsNI)    
    readid<-as.character(reads$id)[queryHits(over)]
    readside<-reads$rid[queryHits(over)]
    exid<-DB@exonsNI$id[subjectHits(over)]
    exst<-start(DB@exonsNI)[subjectHits(over)]
    
    cat("Counting paths\n")
    pCounts<-.Call("pathCounts", readid, readside, exst, exid)
    pCounts<-lapply(pCounts[1:2], function(x) x[1:pCounts[[3]]])
    names(pCounts[[2]])<-pCounts[[1]]
    pCounts<-pCounts[[2]]
    pCounts <- pCounts[grepl("-[0-9]", names(pCounts))]
    
    sel <- strsplit(names(pCounts), split='-|\\.')
    sel1 <- lapply(sel, "[", 2)
    sel1 <- unlist(sel1)
    
    nislEx <- lapply(DB@islands, length)
    nislEx <- rep(names(DB@islands), unlist(nislEx))
    islEx <- unlist(lapply(DB@islands, names))
    names(islEx) <- nislEx
    isl <- match(sel1, islEx)
    isl <- names(islEx)[isl]
    splCounts <- split(pCounts, isl)
    splCounts <- lapply(splCounts, function(x) x[grepl("-", names(x))])
    
    if(DB@denovo){
      sel <- sapply(sel, "[", -1)
      tmp <- split(sel, isl)
      if(mc.cores>1) {
        require(multicore)
        tmp1 <- multicore::mclapply(names(tmp), function(x){
          n <- sapply(tmp[[x]], length)
          nn <- unlist(tmp[[x]])
          names(nn) <- rep(names(splCounts[[x]]), n)
          nnn <- nn %in% names(DB@islands[[x]])
          nnnn <- tapply(nnn, names(nn), all)
          splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
        }, mc.cores=mc.cores)
      } else {
        tmp1 <- lapply(names(tmp), function(x){
          n <- sapply(tmp[[x]], length)
          nn <- unlist(tmp[[x]])
          names(nn) <- rep(names(splCounts[[x]]), n)
          nnn <- nn %in% names(DB@islands[[x]])
          nnnn <- tapply(nnn, names(nn), all)
          splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
        })
      }
      names(tmp1) <- names(tmp)
      splCounts <- tmp1
    }
    
    ans <- vector(length(DB@islands), mode='list')
    names(ans) <- names(DB@islands)
    ans[names(splCounts)] <- splCounts
    ans
  }

genomeBystrand <- function(DB, strand){
  sel <- DB@islandStrand==strand
  islands <- DB@islands[sel]
  transcripts <- DB@transcripts[sel]
  exonsNI <- DB@exonsNI[DB@exonsNI$id %in% unlist(transcripts),]
  exon2island <- DB@exon2island[DB@exon2island$id %in%unlist(transcripts),]
  islandStrand <- DB@islandStrand[sel]
  txid <- unlist(lapply(transcripts, names))
  aliases <- DB@aliases[DB@aliases$tx %in% txid,]
  ans <- new("annotatedGenome", aliases=aliases, denovo=TRUE, exonsNI=exonsNI, islandStrand=islandStrand, transcripts=transcripts, exon2island=exon2island, dateCreated=Sys.Date(), genomeVersion=DB@genomeVersion, islands=islands)
  ans
}

pathCounts<-function(reads, DB, mc.cores=1) {
  if (class(reads) != 'procBam') stop('reads must be an object of class procBam')
  if(!reads@stranded) {
    counts <- procPaths(reads@pbam, DB, mc.cores)
    ans <- new("pathCounts", counts=list(counts), denovo=DB@denovo, stranded=reads@stranded)
  }
  else {
    plusDB <- genomeBystrand(DB, strand="+")
    plus <- procPaths(reads$plus, plusDB, mc.cores)
    minusDB <- genomeBystrand(DB, strand="-")
    minus <- procPaths(reads$minus, minusDB, mc.cores)
    ans <- new("pathCounts", counts=list(plus=plus, minus=minus), denovo=DB@denovo, stranded=reads@stranded)
  }
  ans
}
