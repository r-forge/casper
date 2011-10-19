pathCounts<-function(reads, exons) {
  cat("Finding overlaps between reads and exons\n")
  over<-findOverlaps(reads, exons)
  readid<-as.character(reads$id)[queryHits(over)]
  readside<-reads$rid[queryHits(over)]
  exid<-exons$id[subjectHits(over)]
  exst<-start(exons)[subjectHits(over)]

  cat("Counting paths\n")
  pCounts<-.Call("pathCounts", readid, readside, exst, exid)
  pCounts<-lapply(pCounts[1:2], function(x) x[1:pCounts[[3]]])
  names(pCounts[[2]])<-pCounts[[1]]
  pCounts<-pCounts[[2]]
  pCounts
}
