findNewExons <- function(reads, DB, minReads=1, readLen=NA, stranded=FALSE, pvalFilter=0.05){
  if(is.na(readLen)) stop("No readLen specified")
  if(stranded){
  } else {
#    exons<-RangedData(exons@unlistData)
    exons <- DB$exonsNI
    cat("\tCalculating coverage\n")
    cov<-coverage(reads)
    islands<-slice(cov, lower=1)
    sel<-islands %in% exons
    newisl<-islands[!sel]
    counts<-viewSums(newisl)
    counts<-round(counts/readLen)
    cat("\tGenerating RangedData\n")
    newisl<-RangedData(IRangesList(start=start(newisl), end=end(newisl)), counts=unlist(counts))
    newisl<-RangedData(newisl, counts=jj)
    cat("\tFinding p-values\n")
    n<-sum(newisl$counts)
    p<-1/nrow(newisl)
    pvals<-pbinom(newisl$counts - 1, n, p, lower.tail=FALSE )
    newisl[['pvalue']]<-pvals
    
    newisl <- newisl[newisl[['pvalue']]<=pvalFilter,]
  }
  id<-max(exons$id)+1
  id<-id:(id+nrow(newisl)-1)
  newisl$id<-id
  cat("\tAdding to old exons\n")
  newisl <- RangedData(IRanges(c(start(exons), end(newisl)), c(end(exons), end(newisl))), space=c(as.character(exons$space), as.character(newisl$space)), id=c(exons$id, newisl$id))
  
  return(newisl)
}

makeIslands <- function(allexs){
    
    load("allexs.RData")
    
    exons <- allexs 
    txs <- as.integer(as.factor(names(exons)))
    totEx <- length(exons)
    uniex <- unique(exons);
    nexR <- length(uniex);
    islands <- rep(0, nexR);
    
    dyn.load("casper/pkg/casper/src/casperGamma.so")
    ans<-.Call("makeGeneIslands", exons, islands, uniex, txs, totEx, nexR)
    names(ans) <- uniex
    ans
}

assignExons2Gene <- function(exons, DB, reads, maxDist=5000, stranded=FALSE, mc.cores=1){
    sel<-exons$id %in% DB$exonsNI$id
    newex<-exons$id[!sel]
    
#assign by junctions
    
    cat("Finding overlaps between genes and new exons\n")
    over<-findOverlaps(reads, exons)
    shits <- subjectHits(over)
    qhits<- queryHits(over)
    
#Select new exons
    exs <- exons$id[shits]
    sel <- exs %in% newex
    exs <- exs[sel]
    
    cat("Finding junctions of new exons\n")
#Select read ids in these exons
    rea <- reads$id[qhits]
    rea <- rea[sel]
    
#Select all reads with these ids and overlapping exons
    sel1 <- (1:nrow(reads))[reads$id %in% rea]
    sel2 <- qhits %in% sel1
    rea1 <- reads$id[qhits[sel2]]
    exs1 <- exons$id[shits[sel2]]
    
#Find junctions between exons (new to old and new)
    junx<-tapply(exs1, rea1, unique)
    exids <- rep(1:length(junx), unlist(lapply(junx, length)))
    junx<-unname(unlist(junx))
    names(junx) <- exids
    
#Build gene islands
    txids <- rep(names(DB$newTxs), unlist(lapply(DB$newTxs, length)))
    oldexs <- unlist(DB$newTxs)
    names(oldexs)<-txids
    allexs <- c(oldexs, junx)
    islands <- makeIslands(allexs) 
    
    cat("Formatting denovo genome\nSplitting genes\n")
    exon2gene <- as.data.frame(exons)
    exon2gene$gene <- islands[exon2gene$id]
    rownames(exon2gene) <- exon2gene$id
    
    genes <- split(exon2gene, exon2gene$gene)
    if(mc.cores>1) {
        require(multicore)
        genes<-mclapply(genes, function(x){ y <- IRanges(x$start, x$end); names(y) <- x$id; y}, mc.cores=mc.cores)
    } else {
        genes<-lapply(genes, function(x){ y <- IRanges(x$start, x$end); names(y) <- x$id; y})
    }
    
    cat("Splitting transcripts\n")
    extxs <- unlist(lapply(DB$newTxs, "[", 1))
    sel <- exon2gene[extxs,]$gene
    names(sel) <- names(DB$newTxs)
    transcripts <- DB$newTxs[names(sel)]
    transcripts <- split(transcripts, sel)
    genes <- lapply(names(genes), function(x) list(exons=genes[[x]], transcripts=transcripts[[x]])) 
    
    ans <- list(genes=genes, exon2gene=exon2gene)
    ans
}

      
createDenovoGenome <- function(reads, DB, readLen, stranded, mc.cores=1){
  setClass("denovoGenome", representation(genes = "list", exon2gene= "data.frame"))
  denovo <- new("denovoGenome")
  cat("Finding new exons\n")
  newex <- findNewExons(reads, DB, readLen=readLen, pvalFilter=0.05)
  cat("Done...\nCreating denovo genome\n")
  denovo <- assignExons2Gene(newex, DB, reads, maxDist=5000, stranded=stranded, mc.cores=mc.cores) 
  denovo
}

