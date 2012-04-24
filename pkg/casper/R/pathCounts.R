require(methods)

setClass("pathCounts", representation(counts = "list", denovo = "logical"))

valid_pathCounts <- function(object) {
  msg <- NULL
#validity checks go here
  if(!(is.null(msg))) { TRUE } else { msg }
}

setValidity("pathCounts", valid_pathCounts)
setMethod("show", signature(object="pathCounts"), function(object) {
  if(object@denovo) {
    cat("Denovo pathCounts object with",length(object@counts),"islands and",sum(!unlist(lapply(object@counts, is.null))),"non zero islands.\n")
  } else cat("Known pathCounts object with",length(object@counts),"islands and", sum(!unlist(lapply(object@counts, is.null))),"non zero islands.\n")
})



pathCounts<-function(reads, DB) {

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
  
  sel <- strsplit(names(pCounts), split='-|\\.')
  sel <- lapply(sel, "[", 2)
  sel <- unlist(sel)
  
  nislEx <- lapply(DB@islands, length)
  nislEx <- rep(names(DB@islands), unlist(nislEx))
  islEx <- unlist(lapply(DB@islands, names))
  names(islEx) <- nislEx
  isl <- match(sel, islEx)
  isl <- names(islEx)[isl]

  splCounts <- split(pCounts, isl)

  ans <- vector(length(islands), mode='list')
  names(ans) <- names(islands)
  ans[names(splCounts)] <- splCounts

  ans <- new("pathCounts", counts=ans, denovo=DB@denovo)  
  ans
}
