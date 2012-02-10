calcExp<-function(distrs, genomeDB, pc){
  dexo<-as.data.frame(genomeDB$exonsNI)
  mexo<-as.matrix(dexo[,c(5,4)])
  stafun = ecdf(distrs$stDis)
  fill<-c(rep(0,(min(as.numeric(names(distrs$lenDis))))), distrs$lenDis)
  names(fill)[1:min(as.numeric(names(distrs$lenDis)))]<-0:(min(as.numeric(names(distrs$lenDis)))-1)
  lendis<-as.table(fill/sum(fill))
  exp<-.Call("calc", mexo, genomeDB$newTx, pc, stafun, lendis)
  names(exp)<-names(genomeDB$newTx)
  exp<-as.matrix(exp)
  featureData<-as.numeric(unlist(genomeDB$txs@elementMetadata$gene_id))
  names(featureData)<-unlist(genomeDB$txs@elementMetadata$tx_name)
    featureData<-as.data.frame(as.matrix(featureData[rownames(exp)]))
  metadata <- data.frame(labelDescription = "gene_id", row.names="gene_id")
  featureData<-new("AnnotatedDataFrame", data=featureData, varMetadata=metadata)
  exp<-new("ExpressionSet", exprs=exp, featureData=featureData, annotation="entrez")
  exp
}
