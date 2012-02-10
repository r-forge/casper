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
    

procBam<-function(bamFileName, samtools, chrom, bam, parallel){
  require(IRanges)
  if(!parallel){
    if(bam>0){
      size<-system(paste(samtools, "/samtools idxstats ", bamFileName, sep=""), intern=T)
      size<-do.call(rbind, lapply(strsplit(size, "\t"), function(x) c(x[1], x[3])))
      if(chrom!=""){
        size<-as.numeric(size[size[,1]==chrom,][2])
      } else {
        size<-sum(as.numeric(size[,2]))
      }
    } else {
      if(chrom!=""){      
        size<-system(paste("grep ", chrom, " ", bamFileName, "| wc ", sep=""), intern=T)
      } else size<-system(paste("wc ", bamFileName, sep=""), intern=T)
      size<-strsplit(size, " ")
      size<-as.numeric(size[[1]][size[[1]]!=""][1])
    }
    cat(paste("Reading file with ", size," reads\n"))
    data<-.Call("procBam", bamFileName, samtools,  chrom, bam, size)
    
	data[[1]]<-buildRD(reads=data[[1]], frags=NA)
    data[[2]]<-buildRD(reads=NA, frags=data[[2]])
    names(data)<-c("reads", "frags")
      
  } else {
    if(bam>0){
      size<-system(paste(samtools, "/samtools idxstats ", bamFileName, sep=""), intern=T)
      size<-do.call(rbind, lapply(strsplit(size, "\t"), function(x) c(x[1], x[3])))
      chroms<-size[size[,2]>0,1]
      size<-as.numeric(size[size[,2]>0,2])
      names(size)<-chroms
      data<-lapply(names(size), function(x){
        cat(paste("Reading chromosome ", x," with ", size[x]," reads\n"))
        res<-.Call("procBam", bamFileName, samtools,  x, bam, size[x])
        if(length(res[[1]][[1]])>0){
                   
        res[[1]]<-buildRD(reads=res[[1]], frags=NA)
        res[[2]]<-buildRD(reads=NA, frags=res[[2]])
                   
        names(res)<-c("reads", "frags")
        } else res<-NA
        gc()
        res
      }
      )
      reads<-lapply(data[!is.na(data)], function(x) x$reads)
      frags<-lapply(data[!is.na(data)], function(x) x$frags)
      data<-list()
      data$reads<-do.call(rbind, reads)
      data$frags<-do.call(rbind, frags)
      data
    } else {
      chroms<-system(paste("head -n 100 ", bamFileName, " | grep SQ ", sep=""), intern=T)
      chroms<-lapply(chroms, function(x) strsplit(x, "\t")[[1]][2])
      chroms<-lapply(chroms, function(x) strsplit(x, ":")[[1]][2])
      data<-lapply(chroms, function(x){
			size<-system(paste("awk '{print $3}' ", bamFileName, " | grep ", x, " | wc | awk '{print $1}' ", sep=""), intern=T)
            size<-as.numeric(size)
            cat(paste("Reading chromosome ", x," with ", size," reads\n"))
 	if(size>0){
       res<-.Call("procBam", bamFileName, samtools,  x, bam, size)
	 } else res<-NULL
        if(length(res[[1]][[1]])>0){
            res[[1]]<-buildRD(reads=res[[1]], frags=NA)
                   res[[2]]<-buildRD(reads=NA, frags=res[[2]])
                   names(res)<-c("reads", "frags")

        } else res<-NA
                   gc()
                   res
                   })
        reads<-lapply(data[!is.na(data)], function(x) x$reads)
        frags<-lapply(data[!is.na(data)], function(x) x$frags)
        data<-list()
        data$reads<-do.call(rbind, reads)
        data$frags<-do.call(rbind, frags)
        data
    }
  }
data
}
