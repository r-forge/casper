generateNOexons<-function(exByTx, ExonsU, mc.cores){
    
    exByTx<-unlist(exByTx)
    if(mc.cores>1) require(multicore)
    if(mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        exonsNI<-mclapply(levels(exByTx@seqnames), function(x){
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
    exonsNI$id<-1:nrow(exonsNI)
    overEx<-findOverlaps(exonsNI, ExonsU)
    exid<-exonsNI$id[exonsNI$id[queryHits(overEx)]]
    exkey<-split(exid, exByTx@elementMetadata$exon_id[subjectHits(overEx)])
    if(mc.cores>1) require(multicore)
    if(mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        exkey<-mclapply(exkey, unique, mc.cores=mc.cores)
      } else stop('multicore library has not been loaded!')
    }
    else {
      exkey<-lapply(exkey, unique, mc.cores=mc.cores)
    }
    res<-list(exkey=exkey, exons=exonsNI)
    res 
  }
    

mapEx<-function(exs, exkey){
  nxs<-unname(unlist(exkey[as.character(exs)]))
  nxs
}


procGenome<-function(genome, mc.cores=mc.cores){

  require(GenomicFeatures)
  genDB<-makeTranscriptDbFromUCSC(genome=genome, tablename="refGene")

  cat("Processing Exons and Transcrips\n")
  Exons<-exonsBy(genDB, by="tx")
  txs<-transcripts(genDB,columns=c("tx_id","tx_name","gene_id","exon_id","cds_id"))
  txsRD<-RangedData(txs)

  ExonsDF<-as.data.frame(Exons)
  ExonsRD<-RangedData(ranges=IRanges(start=ExonsDF$start, end=ExonsDF$end), space=ExonsDF$seqnames, element=ExonsDF$element, exon_id=ExonsDF$exon_id, strand=ExonsDF$strand)

  #Find Non-overlapping exons

  cat("Finding non-overlapping exons\n")
  ExonsU<-unlist(Exons)
  fixExons<-generateNOexons(Exons, ExonsU, mc.cores)
  exkey=fixExons$exkey
  exonsNI<-fixExons$exons

  #Find transcript structure for new exons

  cat("Remapping transcript structure to new exons\n")
  
  t<-lapply(txs@elementMetadata$exon_id, unlist)
  s<-as.character(txs@strand)
  exon_ids<-lapply(1:length(t), function(x) if(s[[x]]=="+") t[[x]] else rev(t[[x]]))
  if(mc.cores>1) require(multicore)
  if(mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      newTxs<-multicore::mclapply(exon_ids, function(x) mapEx(x, exkey), mc.cores=mc.cores)
    } else stop('multicore library has not been loaded!')
  } else  {
    newTxs<-lapply(exon_ids, function(x) mapEx(x, exkey))
  }

  names(newTxs)<-txs@elementMetadata$tx_id

  geneids<-txs@elementMetadata$gene_id
  geneids<-unlist(geneids)
  names(geneids)<- names(newTxs)
  txsNewExonID<-list(gene=geneids, exonsNI=exonsNI, exons=Exons, txs=txs, newTxs=newTxs, exonmap=exkey)
  txsNewExonID
 } 
