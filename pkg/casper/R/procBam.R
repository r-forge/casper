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
    
  #process lengths RangedData
    data[[1]]<-lapply(data[[1]][1:length(data[[1]])-1], function(x) x[1:data[[1]][[length(data[[1]])]]])
      if(sum(grepl(unique(data[[1]][[5]]), "chr"))>0) {
          data[[1]]<-RangedData(IRanges(start=ifelse(data[[1]][[1]]<data[[1]][[2]], data[[1]][[1]], data[[1]][[2]]), end=ifelse(data[[1]][[1]]>data[[1]][[2]], data[[1]][[1]], data[[1]][[2]])), space=data[[1]][[5]], id=data[[1]][[4]], flag=data[[1]][[3]], rid=data[[1]][[6]])
      } else {
        data[[1]]<-RangedData(IRanges(start=ifelse(data[[1]][[1]]<data[[1]][[2]], data[[1]][[1]], data[[1]][[2]]), end=ifelse(data[[1]][[1]]>data[[1]][[2]], data[[1]][[1]], data[[1]][[2]])), space=paste("chr", data[[1]][[5]], sep=""), id=data[[1]][[4]], flag=data[[1]][[3]], rid=data[[1]][[6]])
      }
  #Process reads RangedData
      data[[2]]<-lapply(data[[2]][1:length(data[[2]])-1], function(x) x[1:data[[2]][[length(data[[2]])]]])
      data[[2]]<-lapply(data[[2]], function(x) x[data[[2]][[7]]==data[[2]][[8]]])
      if(sum(grepl(unique(data[[2]][[7]]), "chr"))>0) {
    data[[2]]<-RangedData(IRanges(start=ifelse(data[[2]][[1]]<data[[2]][[3]], data[[2]][[1]], data[[2]][[3]]), end=ifelse(data[[2]][[1]]<data[[2]][[3]], data[[2]][[4]], data[[2]][[2]])), data[[2]][[7]], flag1=data[[2]][[5]], flag2=data[[2]][[6]])
      } else {
          data[[2]]<-RangedData(IRanges(start=ifelse(data[[2]][[1]]<data[[2]][[3]], data[[2]][[1]], data[[2]][[3]]), end=ifelse(data[[2]][[1]]<data[[2]][[3]], data[[2]][[4]], data[[2]][[2]])), space=paste("chr", data[[2]][[7]], sep=""), flag1=data[[2]][[5]], flag2=data[[2]][[6]])
      }
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
          res[[1]]<-lapply(res[[1]][1:length(res[[1]])-1], function(x) x[1:res[[1]][[length(res[[1]])]]])
          res[[1]]<-RangedData(IRanges(start=ifelse(res[[1]][[1]]<res[[1]][[2]], res[[1]][[1]], res[[1]][[2]]), end=ifelse(res[[1]][[1]]>res[[1]][[2]], res[[1]][[1]], res[[1]][[2]])), space=paste("chr", res[[1]][[5]], sep=""), id=res[[1]][[4]], flag=res[[1]][[3]], rid=res[[1]][[6]])
  #Process reads RangedData
        res[[2]]<-lapply(res[[2]][1:length(res[[2]])-1], function(x) x[1:res[[2]][[length(res[[2]])]]])
        res[[2]]<-lapply(res[[2]], function(x) x[res[[2]][[7]]==res[[2]][[8]]])
        res[[2]]<-RangedData(IRanges(start=ifelse(res[[2]][[1]]<res[[2]][[3]], res[[2]][[1]], res[[2]][[3]]), end=ifelse(res[[2]][[1]]<res[[2]][[3]], res[[2]][[4]], res[[2]][[2]])), space=paste("chr", res[[2]][[7]], sep=""), flag1=res[[2]][[5]], flag2=res[[2]][[6]])
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
        size<-system(paste("grep ", x, " ", bamFileName, "| wc ", sep=""), intern=T)
        size<-strsplit(size[[1]], "    ")[[1]][2]
        cat(paste("Reading chromosome ", x," with ", size," reads\n"))
        data<-.Call("procBam", bamFileName, samtools,  x, bam, size)
        if(length(data[[1]][[1]])>0){
          data[[1]]<-lapply(data[[1]][1:length(data[[1]])-1], function(x) x[1:data[[1]][[length(data[[1]])]]])
          data[[1]]<-RangedData(IRanges(start=ifelse(data[[1]][[1]]<data[[1]][[2]], data[[1]][[1]], data[[1]][[2]]), end=ifelse(data[[1]][[1]]>data[[1]][[2]], data[[1]][[1]], data[[1]][[2]])), space=paste("chr", data[[1]][[5]], sep=""), id=data[[1]][[4]], flag=data[[1]][[3]], rid=data[[1]][[6]])
          #Process reads RangedData
          data[[2]]<-lapply(data[[2]][1:length(data[[2]])-1], function(x) x[1:data[[2]][[length(data[[2]])]]])
          data[[2]]<-lapply(data[[2]], function(x) x[data[[2]][[7]]==data[[2]][[8]]])
          data[[2]]<-RangedData(IRanges(start=ifelse(data[[2]][[1]]<data[[2]][[3]], data[[2]][[1]], data[[2]][[3]]), end=ifelse(data[[2]][[1]]<data[[2]][[3]], data[[2]][[4]], data[[2]][[2]])), space=paste("chr", data[[2]][[7]], sep=""), flag1=data[[2]][[5]], flag2=data[[2]][[6]])
        } else data<-NA
        reads<-lapply(data[!is.na(data)], function(x) x$reads)
        frags<-lapply(data[!is.na(data)], function(x) x$frags)
        data<-list()
        data$reads<-do.call(rbind, reads)
        data$frags<-do.call(rbind, frags)
        data
      }
      )
    }
  }
data
}
