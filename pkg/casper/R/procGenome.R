require(methods)



## METHODS

setMethod("getIsland",signature(entrezid='character',txid='missing',genomeDB='annotatedGenome'),function(entrezid, txid, genomeDB) {
  as.character(genomeDB@aliases[genomeDB@aliases$gene_id==entrezid,'island_id'][1])
}
)
setMethod("getIsland",signature(entrezid='missing',txid='character',genomeDB='annotatedGenome'),function(entrezid, txid, genomeDB) {
  as.character(genomeDB@aliases[txid,'island_id'])
}
)

setMethod("transcripts", signature(entrezid='missing',islandid='character',genomeDB='annotatedGenome'), function(entrezid, islandid, genomeDB) {
  IRangesList(lapply(genomeDB@transcripts[[islandid]],function(z) genomeDB@islands[[islandid]][as.character(z)]))
}
)
setMethod("transcripts", signature(entrezid='character',islandid='missing',genomeDB='annotatedGenome'), function(entrezid, islandid, genomeDB) {
  txids <- rownames(genomeDB@aliases)[genomeDB@aliases$gene_id==entrezid]
  islandid <- getIsland(entrezid=entrezid,genomeDB=genomeDB)
  IRangesList(lapply(genomeDB@transcripts[[islandid]],function(z) genomeDB@islands[[islandid]][as.character(z)]))[txids]
}
)


### FUNCTIONS

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

generateNOexons<-function(exByTx, startId=1, mc.cores){    
  exByTx<-unlist(exByTx)
    if(mc.cores>1) require(multicore)
    if(mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        exonsNI<-multicore::mclapply(levels(exByTx@seqnames), function(x){
                    y<-exByTx[exByTx@seqnames==x,]
                    st<-start(y)
                    en<-end(y)
                    both<-sort(unique(c(st,en)))
                    rd<-IRanges(start=both[-length(both)], end=both[-1])
                    strd<-start(rd)
                    enrd<-end(rd)
                    strd[strd %in% en]<-strd[strd %in% en]+1
                    enrd[enrd %in% st]<-enrd[enrd %in% st]-1
                    allRD<-IRanges(start=strd, end=enrd)
                  }, mc.cores=mc.cores
                 )
      } else stop('multicore library has not been loaded!')
    }
    else {
      exonsNI<-lapply(levels(exByTx@seqnames), function(x){
        y<-exByTx[exByTx@seqnames==x,]
        st<-start(y)
        en<-end(y)
        both<-sort(unique(c(st,en)))
        rd<-IRanges(start=both[-length(both)], end=both[-1])
        strd<-start(rd)
        enrd<-end(rd)
        strd[strd %in% en]<-strd[strd %in% en]+1
        enrd[enrd %in% st]<-enrd[enrd %in% st]-1
        allRD<-IRanges(start=strd, end=enrd)
      })
       }
    names(exonsNI)<-levels(exByTx@seqnames)
    exonsNI<-IRangesList(exonsNI)
    exonsNI<-RangedData(exonsNI)
    exonsNI <- subsetByOverlaps(exonsNI, exByTx) 
    exonsNI$id<-startId:(nrow(exonsNI)+startId-1)
    overEx<-findOverlaps(exonsNI, exByTx)
    exid<-exonsNI$id[queryHits(overEx)]
    exkey<-split(exid, values(exByTx)$exon_id[subjectHits(overEx)])
    if(mc.cores>1) require(multicore)
    if(mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        exkey<-multicore::mclapply(exkey, unique, mc.cores=mc.cores)
      } else stop('multicore library has not been loaded!')
    }
    else {
      exkey<-lapply(exkey, unique, mc.cores=mc.cores)
    }
    res<-list(exkey=exkey, exons=exonsNI)
    res 
  }

mapEx<-function(exs, exkey){
  nxs<-unname(unlist(exkey[as.character(sprintf("%d", exs))]))
  nxs
}


procGenome<-function(genome, mc.cores=1){

  genDB<-makeTranscriptDbFromUCSC(genome=genome, tablename="refGene")
  cat("Processing Exons and Transcrips\n")
  txs<-transcripts(genDB,columns=c("tx_id","tx_name","gene_id","exon_id","cds_id"))
  txs<-txs[match(unique(unlist(txs@elementMetadata$tx_name)), unlist(txs@elementMetadata$tx_name)),]
  exid <- sapply(txs@elementMetadata$exon_id, function(x) paste(unlist(x), collapse="."))
  names(exid) <- values(txs)$tx_name
  
  #aliases <- txs[exid %in% names(table(exid))[table(exid)>1],]
  txs <- txs[match(unique(exid), exid),]
  
  aliases <- values(txs)[1:3]
  aliases[,3] <- unlist(aliases[,3])
  aliases <- as.data.frame(aliases, stringsAsFactors=F)
  aliases <- cbind(aliases, exid=exid[as.character(aliases$tx_name)])
  alitx <- split(as.character(aliases[,2]), aliases[,4])
  names(alitx) <- unlist(lapply(alitx, function(x) ifelse(length(x)>1, x[x %in% txs@elementMetadata$tx_name], x)))
  nalitx <- rep(names(alitx), unlist(lapply(alitx, length)))
  alitx <- unlist(alitx)
  names(alitx) <- nalitx
  aliases <- cbind(aliases, tx=names(alitx)[match(aliases[,2], alitx)])

  Exons<-exonsBy(genDB, by="tx")
  Exons<-Exons[names(Exons) %in% txs@elementMetadata$tx_id,]
  txnames<-match(names(Exons), txs@elementMetadata$tx_id)
  names(Exons)<-txs@elementMetadata$tx_name[txnames]
  txsRD<-RangedData(txs)
  ExonsDF<-as.data.frame(Exons)
  ExonsRD<-RangedData(ranges=IRanges(start=ExonsDF$start, end=ExonsDF$end), space=ExonsDF$seqnames, element=ExonsDF$element, exon_id=ExonsDF$exon_id, strand=ExonsDF$strand)
  #Find Non-overlapping exons

  cat("Finding non-overlapping exons\n")

  txStrand <- as.character(txs@strand)
  names(txStrand) <- values(txs)$tx_name
  sel <- names(txStrand)[txStrand=="+"]; plusExons <- Exons[sel]
  sel <- names(txStrand)[txStrand=="-"]; minusExons <- Exons[sel]

  fixExonsP <- generateNOexons(plusExons, startId=1, mc.cores)
  fixExonsM <- generateNOexons(minusExons, startId=max(fixExonsP$exons$id)+1, mc.cores)
  
  exkey <- fixExonsP$exkey
  exkey <- c(exkey, fixExonsM$exkey)
  exonsNI <- as.data.frame(fixExonsP$exons)
  exonsNI <- rbind(exonsNI, as.data.frame(fixExonsM$exons))
  exonsNI <- RangedData(exonsNI)
  
  #Find transcript structure for new exons

  cat("Remapping transcript structure to new exons\n")
  s<-ifelse(as.character(Exons@unlistData@strand)=="+", 1, -1)
  exids <- values(Exons@unlistData)$exon_id
  idx <- values(Exons@unlistData)$exon_rank
  idx <- cumsum(idx==1)
  exids <- s*exids
  exon_ids <- tapply(exids, idx, function(x) if(any(x<0)) {rev(-x)} else x)
  names(exon_ids) <- names(Exons)
  
  
  if(mc.cores>1) require(multicore)
  if(mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      newTxs<-multicore::mclapply(exon_ids, function(x) mapEx(x, exkey), mc.cores=mc.cores)
    } else stop('multicore library has not been loaded!')
  } else  {
    newTxs<-lapply(exon_ids, function(x) mapEx(x, exkey))
  }
  names(newTxs)<-names(Exons)

  geneids<-txs@elementMetadata$gene_id
  geneids<-unlist(geneids)
  names(geneids)<- unlist(txs@elementMetadata$tx_id)
  # Make islands

  ex2tx <- unlist(newTxs)
  names(ex2tx) <- rep(names(newTxs), unlist(lapply(newTxs, length)))
  islands <- makeIslands(ex2tx)
  
  exon2island <- as.data.frame(exonsNI)
  exon2island$island <- islands[as.character(exon2island$id)]
  rownames(exon2island) <- exon2island$id
  islands <- split(exon2island, exon2island$island)
  txStrand <- as.character(strand(txs))
  names(txStrand) <- values(txs)$tx_name
  extxs <- unlist(lapply(newTxs, "[", 1))  

  
  if(mc.cores>1) require(multicore)
  if(mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      islands <- multicore::mclapply(islands, function(x){
        y <- IRanges(x$start, x$end);
        names(y) <- x$id;
        if(txStrand[names(ex2tx)[ex2tx %in% x$id[1]]]=="-") y <- rev(y)
        y
      }, mc.cores=mc.cores) } else stop('multicore library has not been loaded!')
  } else {
    islands <- lapply(islands, function(x){ y <- IRanges(x$start, x$end); names(y) <- x$id; y})
  }

  cat("Splitting transcripts\n")

  sel <- match(extxs, exon2island$id)
  tx2island <- exon2island$island[sel]
  names(tx2island) <- names(newTxs)
  transcripts <- newTxs[names(tx2island)]
  transcripts <- split(transcripts, tx2island)
  sel <- unlist(lapply(transcripts, function(x) names(x)[1]))
  islandStrand <- txStrand[sel]
  names(islandStrand) <- names(islands)

  sel <- islandStrand=='-'
  transcripts[sel] <- lapply(transcripts[sel], function(x) lapply(x,rev))

  id2tx <- data.frame(island_id=rep(names(transcripts),sapply(transcripts,length)) , txname=unlist(sapply(transcripts,names)))
  rownames(id2tx) <- id2tx$txname
  aliases$island_id <- id2tx[rownames(aliases),'island_id']

  ans <- new("annotatedGenome", islands=islands, transcripts=transcripts, exon2island=exon2island, aliases=aliases, exonsNI=exonsNI, islandStrand=islandStrand, dateCreated=Sys.Date(), genomeVersion=genome, denovo=FALSE)
  ans
 } 

