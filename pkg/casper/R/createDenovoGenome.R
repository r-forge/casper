require(methods)

findNewExons <- function(reads, DB, minReads=1, readLen=NA, stranded=FALSE, pvalFilter=0.05){
  if(is.na(readLen)) stop("No readLen specified")
  if(stranded){
  } else {
    exons <- DB@exonsNI
    #cat("\tCalculating coverage\n")
    cov<-coverage(reads)
    islands<-slice(cov, lower=1)
    sel<-islands %in% exons
    if(all(sel)) {
      cat("No new islands found\n")
      return(NULL)
    }
    newisl<-islands[!sel]
    counts<-viewSums(newisl)
    counts<-round(counts/readLen)
    #cat("\tGenerating RangedData\n")
    newisl<-RangedData(IRangesList(start=start(newisl), end=end(newisl)), counts=unlist(counts))
    #cat("\tFinding p-values\n")
    n<-sum(newisl$counts)
    p<-1/nrow(newisl)
    pvals<-pbinom(newisl$counts - 1, n, p, lower.tail=FALSE )
    newisl[['pvalue']]<-pvals
    newisl <- newisl[newisl[['pvalue']]<=pvalFilter,]
  }
  id<-max(exons$id)+1
  id<-id:(id+nrow(newisl)-1)
  newisl$id<-id
  #cat("\tAdding to old exons\n")
  newisl <- RangedData(IRanges(c(start(exons), start(newisl)), c(end(exons), end(newisl))), space=c(as.character(exons$space), as.character(newisl$space)), id=c(exons$id, newisl$id))
  sel <- !(newisl %in% DB@exonsNI)
  newisl <- newisl[sel,]
  return(newisl)
}

makeIslands <- function(allexs){
    exons <- allexs
      txs <- as.integer(as.factor(names(exons)))
      totEx <- length(exons)
      uniex <- unique(exons)
      nexR <- length(uniex)
      islands <- rep(0, nexR)
      ans<-.Call("makeGeneIslands", exons, islands, uniex, txs, totEx, nexR)
      names(ans) <- uniex
      ans
  }

assignExons2Gene <- function(exons, DB, reads, maxDist=1000, minLinks=2, maxLinkDist=100000, stranded=FALSE, mc.cores=1){

  newex <- exons$id
  exons <- RangedData(IRanges(start=c(start(DB@exonsNI), start(exons)), end=c(end(DB@exonsNI), end(exons))), space=c(as.character(space(DB@exonsNI)), as.character(space(exons))), id=c(DB@exonsNI$id, exons$id))
  over <- findOverlaps(reads, exons)
  shits <- subjectHits(over)
  qhits <- queryHits(over)
  
#Select new exons
  exs <- exons$id[shits]
  sel <- exs %in% newex
  exs <- exs[sel]
#Select read ids in these exons
  rea <- reads$id[qhits]
  rea <- rea[sel]

  sel1 <- (1:nrow(reads))[reads$id %in% rea]
  sel2 <- qhits %in% sel1
  rea1 <- reads$id[qhits[sel2]]
  exs1 <- exons$id[shits[sel2]]
  len <- length(unique(rea1))
  if(length(unique(exs1))>1)  ans <- .Call("joinExons", exs1, rea1, len)
  junx <- ans[[1]][ans[[2]]>=minLinks]
  junx <- strsplit(junx, split=".", fixed=T)
  names(junx) <- 1:length(junx)
  nalljunx <- rep(1:length(junx), unlist(lapply(junx, length)))
  alljunx <- unlist(junx)
  names(alljunx) <- nalljunx

  tmpex <- as.data.frame(exons)
  rownames(tmpex) <- tmpex$id
  tmpex <- tmpex[alljunx,]
  pos <- (tmpex$end+tmpex$start)/2
  dis <- tapply(pos, names(alljunx), function(x) (max(x)-min(x) < maxLinkDist ))
  alljunx <- alljunx[names(alljunx) %in% as.character(names(dis)[dis])]
  sing <- newex[!(newex %in% unique(alljunx))]
  mname <- max(as.numeric(names(alljunx)))
  names(sing) <- (mname+1):(mname+length(sing)) 
  alljunx <- c(alljunx, sing)

  cat("Building islands\n")
#Build gene islands
  oldexs <- unlist(DB@transcripts, recursive=F)
  txids <- sub("[0-9]+\\.", "", names(oldexs))
  txids <- rep(txids, unlist(lapply(oldexs, length)))
  oldexs <- unlist(oldexs)
  names(oldexs) <- txids
  allexs <- c(oldexs, alljunx)
  islands <- makeIslands(allexs) 
  exons$island <- islands[as.character(exons$id)] 
  exon2gene <- as.data.frame(exons)
  exon2gene$gene <- islands[as.character(exon2gene$id)]
  rownames(exon2gene) <- exon2gene$id
  islands <- split(exon2gene, exon2gene$gene)

  if(mc.cores>1) require(multicore)
  if(mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      islands<-multicore::mclapply(islands, function(x){ y <- IRanges(x$start, x$end); names(y) <- x$id; y}, mc.cores=mc.cores)
    } else stop('multicore library has not been loaded!')
  }  else {
    islands<-lapply(islands, function(x){ y <- IRanges(x$start, x$end); names(y) <- x$id; y})
  }

#Link islands by distance

  if(maxDist>0){

    cat("\tJoining islands by distance\n")
    old <- exons[!(exons$id %in% newex),]
    new <- exons[exons$id %in% newex,]
    
    mapisl <- lapply(names(exons), function(i){
      if(sum(space(new)==i)>0 & sum(space(old)==i)>0){

        isl <- by(as.data.frame(ranges(new)[[i]]), new[i]$island, function(x) c(min(x$start), max(x$end)))
        nisl <- names(isl)
        isl <- do.call(rbind, isl)
        isl <- IRanges(isl[,1], isl[,2])
        ne <- distanceToNearest(isl, ranges(old)[[i]])
        ne$isl <- nisl[ne$query]
        sel <- ne$distance < maxDist
        ne <- ne[sel,]
        ne$oisl <- old[i]$island[ne$subject]
        ne <- ne[ne$isl != ne$oisl,]
        rownames(ne) <- new[i]$id[ne$query]
      } else ne <- NA
      ne
    })
    names(mapisl)<-names(exons)
    
    newexons <- lapply(names(new), function(i){
      tmp <- new[i]
      if(length(mapisl[[i]])>1){
        sel <- as.character(tmp$id) %in% rownames(mapisl[[i]])
        map <- match(as.character(tmp$id)[sel], rownames(mapisl[[i]]))
        tmp$island[sel] <- mapisl[[i]]$oisl[map]      
      }
    tmp
    }
                       )
    names(newexons) <- names(new)
  
    exon2island <- rbind(as.data.frame(old), newexons)
    rownames(exon2island) <- exon2island$id
  }

  exon2island <- as.data.frame(exons)
  
  extxs <- unlist(DB@transcripts, recursive=F)
  names(extxs) <- sub("[0-9]+\\.", "", names(extxs))
  sel <- unlist(lapply(extxs, "[", 1))
  sel <- match(sel, exon2island$id)
  tx2gene <- exon2gene$island[sel]
  transcripts <- vector(length=length(islands), mode="list")
  names(transcripts) <- names(islands)
  tmp <- split(extxs, tx2gene)
  transcripts[names(tmp)] <- tmp

  islands <- split(exon2island, exon2island$island)
  cat("fixing genome structure\n")
  exStrand <- rep(DB@islandStrand, unlist(lapply(DB@islands, length)))
  names(exStrand) <- unlist(lapply(DB@islands, names))
  sel <- unlist(lapply(islands, function(x) x$id[1]))
  names(sel) <- names(islands)
  sel2 <- exStrand[as.character(sel)]
  names(sel2) <- names(sel)
  sel2 <- sel2[!is.na(sel2)]
  islandStrand <- vector(length=length(islands))
  names(islandStrand) <- names(islands)
  islandStrand[names(sel2)] <- sel2
  sel <- islandStrand=="-"
  islandStrand[islandStrand=="FALSE"] <- NA
  
  if(mc.cores>1) {
    require(multicore)
    if(unique(islandStrand[!is.na(islandStrand)])=='+') {
      islands <- multicore::mclapply(names(islands), function(x) {ord <- order(islands[[x]]$start); y <- IRanges(islands[[x]]$start[ord], islands[[x]]$end[ord]); names(y) <- islands[[x]]$id[ord]; y}, mc.cores=mc.cores)
    } else islands <- multicore::mclapply(names(islands), function(x) {ord <- order(islands[[x]]$start, decreasing=T); y <- IRanges(islands[[x]]$start[ord], islands[[x]]$end[ord]); names(y) <- islands[[x]]$id[ord]; y}, mc.cores=mc.cores)
  } else {
    if(unique(islandStrand[!is.na(islandStrand)])=='+') {
      islands <- lapply(names(islands), function(x) {ord <- order(islands[[x]]$start); y <- IRanges(islands[[x]]$start[ord], islands[[x]]$end[ord]); names(y) <- islands[[x]]$id[ord]; y})
    } else islands <- lapply(names(islands), function(x) {ord <- order(islands[[x]]$start, decreasing=T); y <- IRanges(islands[[x]]$start[ord], islands[[x]]$end[ord]); names(y) <- islands[[x]]$id[ord]; y}, mc.cores=mc.cores)
  }

  ans <- new("annotatedGenome", islands=islands, transcripts=transcripts, exon2island=exon2island, exonsNI=exons, islandStrand=islandStrand, aliases=DB@aliases, genomeVersion=DB@genomeVersion, dateCreated=Sys.Date(), denovo=TRUE)
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

mergeStrDenovo <- function(plus, minus){  
  nullplus <- unlist(lapply(plus@transcripts, is.null))
  newplus <- unlist(lapply(plus@islands[nullplus], function(x) paste(names(x), collapse=".")))
  #names(newplus) <- names(plus@islands)[nullplus]
  names(minus@islands) <- 1:length(minus@islands)
  nullminus <- unlist(lapply(minus@transcripts, is.null))
  newminus <- unlist(lapply(minus@islands[nullminus], function(x) paste(sort(names(x)), collapse=".")))
  names(newminus) <- names(minus@islands)[nullminus]
  common <- names(newminus)[!(newminus %in% newplus)]
  #common <- newminus %in% newplus
  allislands <- c(plus@islands, minus@islands[!(names(minus@islands) %in% common)])
  names(allislands) <- 1:length(allislands)
  allstrand <- c(plus@islandStrand, minus@islandStrand[!(names(minus@islands) %in% common)])
  names(allstrand) <- 1:length(allstrand)
  alltrans <- c(plus@transcripts, minus@transcripts[!(names(minus@islands) %in% common)])
  names(alltrans) <- 1:length(alltrans)
  allexonsNI <- rbind(as.data.frame(plus@exonsNI), as.data.frame(minus@exonsNI))
  allexonsNI <- allexonsNI[match(unique(allexonsNI$id), allexonsNI$id),]
  ex2is <- allexonsNI
  islen <- lapply(allislands, length)
  islnam <- rep(names(allislands), islen)
  names(islnam) <- unlist(lapply(allislands, names))
  ex2is$island <- islnam[as.character(ex2is$id)]
  allexonsNI <- RangedData(allexonsNI)
  ans <- new("annotatedGenome", aliases=plus@aliases, denovo=TRUE, exonsNI=allexonsNI, islandStrand=allstrand, transcripts=alltrans, exon2island=ex2is, dateCreated=Sys.Date(), genomeVersion=plus@genomeVersion, islands=allislands)
  ans
  }

createDenovoGenome <- function(reads, DB, readLen, stranded=FALSE,  minLinks=2, maxLinkDist=1e+05, maxDist=1000, mc.cores=1){
  if (missing(readLen)) stop('readLen must be specified')
  cat("Finding new exons\n")
  somex <- NULL
  if(reads$stranded){
    newexplus <- findNewExons(reads$plus, DB, readLen=readLen, pvalFilter=0.05)
    newexminus <- findNewExons(reads$minus, DB, readLen=readLen, pvalFilter=0.05)
    if(!is.null(newexplus) | !is.null(newexminus)) somex <- 1
  } else {
    newex <- findNewExons(reads$pbam, DB, readLen=readLen, pvalFilter=0.05)
    if(!is.null(newex)) somex <- 1
  }
  
  if(!is.null(somex)){
    #reads$id <- cumsum(reads$id == c(reads$id[-1], reads$id[1]))
    cat("Done...\nCreating denovo genome for positive strand\n")
    DBplus <- genomeBystrand(DB, "+")
    if(reads$stranded) denovoplus <- assignExons2Gene(newexplus, DBplus, reads$plus, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
    else denovoplus <- assignExons2Gene(newex, DBplus, reads$pbam, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
    cat("Done...\nCreating denovo genome for negative strand\n")
    DBminus <- genomeBystrand(DB, "-")
    if(reads$stranded) denovominus <- assignExons2Gene(newexminus, DBminus, reads$plus, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
    else denovominus <- assignExons2Gene(newex, DBminus, reads$pbam, maxDist=maxDist, stranded=stranded, minLinks=minLinks, maxLinkDist=maxLinkDist, mc.cores=mc.cores)
    cat("Done...\nMerging denovo genome\n")
    denovo <- mergeStrDenovo(denovoplus, denovominus)
  } else {
    denovo <- DB
    denovo@denovo <- TRUE
  }
  denovo

}

 
