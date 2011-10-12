pathCounts<-function(reads, exons) {
  cat("Finding overlaps between reads and exons\n")
  over<-findOverlaps(reads, exons)
  cat("Preparing data\n")
  reads<-reads[queryHits(over),]
  readid<-as.character(reads$id)
  readside<-reads$rid
  exons<-exons[subjectHits(over),]
  exid<-exons$id
  exst<-start(exons)

  cat("Counting paths\n")
  pCounts<-.Call("pathCounts", readid, readside, exst, exid)
  pCounts<-lapply(pCounts[1:2], function(x) x[1:pCounts[[3]]])
  names(pCounts[[2]])<-pCounts[[1]]
  pCounts<-pCounts[[2]]
  pCounts
}
