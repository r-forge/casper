procBam<-function(bamFileName, samtools, chrom){
  require(IRanges)
  data<-.Call("procBam", bamFileName, samtools,  chrom)

  #process lengths RangedData
  data[[1]]<-lapply(data[[1]][1:length(data[[1]])-1], function(x) x[1:data[[1]][[length(data[[1]])]]])
  data[[1]]<-RangedData(IRanges(start=ifelse(data[[1]][[1]]<data[[1]][[2]], data[[1]][[1]], data[[1]][[2]]), end=ifelse(data[[1]][[1]]>data[[1]][[2]], data[[1]][[1]], data[[1]][[2]])), space=paste("chr", data[[1]][[5]], sep=""), id=data[[1]][[4]], flag=data[[1]][[3]], rid=data[[1]][[6]])
  
  #Process reads RangedData
  data[[2]]<-lapply(data[[2]][1:length(data[[2]])-1], function(x) x[1:data[[2]][[length(data[[2]])]]])
  data[[2]]<-lapply(data[[2]], function(x) x[data[[2]][[7]]==data[[2]][[8]]])
  data[[2]]<-RangedData(IRanges(start=ifelse(data[[2]][[1]]<data[[2]][[3]], data[[2]][[1]], data[[2]][[3]]), end=ifelse(data[[2]][[1]]<data[[2]][[3]], data[[2]][[4]], data[[2]][[2]])), space=paste("chr", data[[2]][[7]], sep=""), flag1=data[[2]][[5]], flag2=data[[2]][[6]])
  names(data)<-c("reads", "frags")
  data
}
