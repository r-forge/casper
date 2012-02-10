buildRD<-function(reads, frags){
    if(length(reads)>1){
        reads<-lapply(reads[1:length(reads)-1], function(x) x[1:reads[[length(reads)]]])
        if(sum(grepl("chr", unique(reads[[5]])))>0) {
            reads<-RangedData(IRanges(start=ifelse(reads[[1]]<reads[[2]], reads[[1]], reads[[2]]), end=ifelse(reads[[1]]>reads[[2]], reads[[1]], reads[[2]])), space=reads[[5]], id=reads[[4]], flag=reads[[3]], rid=reads[[6]])
        } else {
            reads<-RangedData(IRanges(start=ifelse(reads[[1]]<reads[[2]], reads[[1]], reads[[2]]), end=ifelse(reads[[1]]>reads[[2]], reads[[1]], reads[[2]])), space=paste("chr", reads[[5]], sep=""), id=reads[[4]], flag=reads[[3]], rid=reads[[6]])
        }
        res<-reads
    } else {
        frags<-lapply(frags[1:length(frags)-1], function(x) x[1:frags[[length(frags)]]])
        frags<-lapply(frags, function(x) x[frags[[7]]==frags[[8]]])
        if(sum(grepl("chr", unique(frags[[7]])))>0) {
            frags<-RangedData(IRanges(start=ifelse(frags[[1]]<frags[[3]], frags[[1]], frags[[3]]), end=ifelse(frags[[1]]<frags[[3]], frags[[4]], frags[[2]])), space=frags[[7]], flag1=frags[[5]], flag2=frags[[6]])
        } else {
            frags<-RangedData(IRanges(start=ifelse(frags[[1]]<frags[[3]], frags[[1]], frags[[3]]), end=ifelse(frags[[1]]<frags[[3]], frags[[4]], frags[[2]])), space=paste("chr", frags[[7]], sep=""), flag1=frags[[5]], flag2=frags[[6]])
        }
        res<-frags
    }
    res
}
    
nbReads <- function(bam0) {
    tab <- table(bam0$cigar)
    count <- sapply(gregexpr("M",names(tab)),length)
    sum(tab*count)
}

procBam<-function(bam){
    require(IRanges)
 
    cat("Calculating total number of reads...\n")
    nreads<-nbReads(bam)
    
    cat("done.\nProcessing cigars and building read's object...\n")
    data<-.Call("procBam", bam$qname, bam$flag,  bam$rname, bam$pos, bam$cigar, length(bam$pos), nreads)
    cat("done.\n")
	data[[1]]<-buildRD(reads=data[[1]], frags=NA)
    data[[2]]<-buildRD(reads=NA, frags=data[[2]])
    names(data)<-c("reads", "frags")
    data
}
