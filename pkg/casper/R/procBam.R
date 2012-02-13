buildRD<-function(reads){
    if(sum(grepl("chr", unique(reads[[5]])))>0) {
        reads<-RangedData(IRanges(start=ifelse(reads[[1]]<reads[[2]], reads[[1]], reads[[2]]), end=ifelse(reads[[1]]>reads[[2]], reads[[1]], reads[[2]])), space=reads[[5]], id=reads[[4]], flag=reads[[3]], rid=reads[[6]])
    } else {
        reads<-RangedData(IRanges(start=ifelse(reads[[1]]<reads[[2]], reads[[1]], reads[[2]]), end=ifelse(reads[[1]]>reads[[2]], reads[[1]], reads[[2]])), space=paste("chr", reads[[5]], sep=""), id=reads[[4]], flag=reads[[3]], rid=reads[[6]])
    }
    res<-reads
    res
}

uniquifyQname<-function(bam, seed){
    sel<-bam$pos>bam$mpos
    bam$qname[sel]<-paste(bam$qname[sel], ".", bam$pos[sel], ".",bam$mpos[sel], sep="")
    sel<-bam$pos<bam$mpos
    bam$qname[sel]<-paste(bam$qname[sel], ".", bam$mpos[sel], ".",bam$pos[sel], sep="")
# Randomly choose "multihits" with same start position but different cigar                                                                               
    tab<-table(bam$qname)
    probPos<-bam$qname %in% names(tab)[tab>2]
    probID<-paste(bam$qname[probPos], bam$pos[probPos])
    set.seed(seed)
    sel<-unlist(tapply(1:sum(probPos), probID, function(x) x[sample(1:length(x), 1)]))
    fixed<-lapply(lapply(bam, "[", probPos), "[", sel)
    res<-lapply(names(bam), function(x) c(bam[[x]], fixed[[x]]))
    names(res)<-names(bam)
    res
}

    
nbReads <- function(bam0) {
    tab <- table(bam0$cigar)
    count <- sapply(gregexpr("M",names(tab)),length)
    sum(tab*count)
}

procBam<-function(bam, seed=1){
    require(IRanges)
    bam<-uniquifyQname(bam, seed)
    cat("Calculating total number of reads...\n")
    nreads<-nbReads(bam)    
    cat("done.\nProcessing cigars and building read's object...\n")
    data<-.Call("procBam", bam$qname, bam$flag,  bam$rname, bam$pos, bam$cigar, length(bam$pos), nreads)
    cat("done.\n")
	ans<-buildRD(data[[1]])
    ans
}
