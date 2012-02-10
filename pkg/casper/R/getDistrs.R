startsReadsTx<-function(ovfr, frags, tx, exonsRD){
  pos<-exonsRD$exon_id %in% unlist(tx$exon_id)
  ss<-start(exonsRD)[pos]
  ee<-end(exonsRD)[pos]
  len<-sum(width(exonsRD)[pos])
  str<-as.character(tx$strand)
  if(str=="-") {
    ss<-rev(ss)
    ee<-rev(ee)
  }
  if(!is.unsorted(ss)){
    pos<-unlist(ovfr[as.character(tx$exon_id)])
    stff<-start(frags)[pos]
    if(length(stff)>0){
      or<-order(stff)
      stff<-stff[stff>=min(ss)]
      stff<-stff[stff<=max(ee)]
      gaps<-ss[-1]-ee[-length(ss)]
      gaps<-c(0,gaps)
      qq<-findInterval(stff, ss)
      res<-((stff-cumsum(gaps)[qq])-min(ss))/len
      wi=width(frags)[pos][or]	
      if(str=="-") {
        res<-rev(res)
        wi<-rev(wi)
      }
      res<-list(res=res, wi=wi)
    } else res=NA
  } else res=NA
  res
}

getDistrs<-function(txs, exons, frags, mc.cores=1){
  
  #Find fragment length distribution for fragments aligning to exons larger than 1000 bases
  cat("Calculating fragment length distribution\n")
  exonsRD<-as.data.frame(exons@unlistData)
  exonsRD<-exonsRD[,1:6]
  colnames(exonsRD)[1]<-"space"
  exonsRD<-RangedData(exonsRD)
  fragsL<-frags[space(frags) %in% space(exonsRD),]
  fragsL<-subsetByOverlaps(fragsL, exonsRD[width(exonsRD)>1000,], type="within")
  ld<-width(fragsL)
  ld<-table(ld)

  #Find fragment start distribution for fragments aligning to transcripts in genes with only one annotated tx

  cat("Calculating fragment start distribution\n")
  txs<-RangedData(txs)
  oneTx<- unlist(txs$gene_id) %in% names(which(table(unlist(txs$gene_id))==1))
  oneTx<-txs[oneTx,]
  oneTx<-subsetByOverlaps(oneTx, frags)
  exonsRD<-subsetByOverlaps(exonsRD, oneTx)
  over<-findOverlaps(frags, exonsRD, type="within")
  ovfr<-split(queryHits(over), exonsRD$exon_id[subjectHits(over)])
  
  if(mc.cores>1) require(multicore)
  if(mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      stPerTx<-mclapply(1:nrow(oneTx), function(x) startsReadsTx(ovfr, frags, oneTx[x,], exonsRD) , mc.cores=mc.cores)
    } else stop('multicore library has not been loaded!')
  } else  {
    stPerTx<-lapply(1:nrow(oneTx), function(x) startsReadsTx(ovfr, frags, oneTx[x,], exonsRD))
  }

  stPerTx<-stPerTx[!is.na(stPerTx)]
  res<-lapply(stPerTx, function(x) x$res)
  wi<-lapply(stPerTx, function(x) x$wi)
  
  res<-unlist(res)
  list(lenDis=ld, stDis=res)
}
