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

pathCounts<-function(reads, DB, mc.cores=1) {

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
  pCounts <- pCounts[grepl("-", names(pCounts))]
  
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
      tmp1 <- mclapply(names(tmp), function(x){
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

  ans <- new("pathCounts", counts=ans, denovo=DB@denovo)  
  ans
}
