buildGene<-function(genomeDB, gene){
    rd<-genomeDB$exons[genomeDB$txs[unlist(genomeDB$txs@elementMetadata$gene_id) == gene]@elementMetadata$tx_name]
    rdl<-lapply(rd, ranges)
    rdl<-IRangesList(rdl)
    rdl
}

genPlot<-function(goi, genomeDB, reads, exp){
  txexp<-exprs(exp)[as.character(genomeDB$txs[unlist(genomeDB$txs@elementMetadata$gene_id) == goi]@elementMetadata$tx_name),]
  if(sum(txexp)>0) txexp<-(txexp/sum(txexp))*100
  gene<-buildGene(gene=goi, genomeDB)
  rangesPlot(x=reads[unique(as.character(seqnames(genomeDB$txs[unlist(genomeDB$txs@elementMetadata$gene_id) == goi])))], gene)
  gene<-RangedData(unlist(gene), space=as.character(unique(seqnames(genomeDB$txs[unlist(genomeDB$txs@elementMetadata$gene_id) == goi]))))
  list(gene=gene, exp=txexp)
    
}

setMethod("genePlot",signature(gene='IRanges'),
  function(gene, genome, refflat, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', ...) {
    if (missing(xlim)) xlim <- range(c(start(gene),end(gene)))
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- 1
    y0 <- seq(1/(nsplice+1),1-1/(1+nsplice),by=1/(1+nsplice))
    x0 <- start(gene)
    x1 <- end(gene)
    if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y=y0,y1=y0)
    segments(x0=x0,x1=x1,y0=y0+.33/(1+nsplice),y1=y0+.33/(1+nsplice))
    segments(x0=x0,x1=x1,y0=y0-.33/(1+nsplice),y1=y0-.33/(1+nsplice))
    segments(x0=x0,x1=x0,y0=y0-.33/(1+nsplice),y1=y0+.33/(1+nsplice))
    segments(x0=x1,x1=x1,y0=y0-.33/(1+nsplice),y1=y0+.33/(1+nsplice))
  }
)

setMethod("genePlot",signature(gene='RangedData'),
  function(gene, genome, refflat, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', ...) {
    genePlot(ranges(gene), genome=genome, refflat=refflat, names.arg=names.arg, xlab=xlab, ylab=ylab, xlim=xlim, cex=cex, yaxt=yaxt, ...)
  }
)

setMethod("genePlot",signature(gene='IRangesList'),
  function(gene, genome, refflat, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', ...) {
    if (missing(xlim)) xlim <- range(c(start(gene),end(gene)))
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- length(gene)
    cols<-rep(rainbow(min(nsplice, 10)), ceiling(nsplice/10))
    y0 <- seq(1/(nsplice+1),1-1/(1+nsplice),by=1/(1+nsplice))
    for (i in 1:nsplice) {
      x0 <- start(gene[[i]])
      x1 <- end(gene[[i]])
      if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y=y0[i],y1=y0[i], col=cols[i])
      segments(x0=x0,x1=x1,y0=y0[i]+.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=cols[i])
      segments(x0=x0,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]-.33/(1+nsplice), col=cols[i])
      segments(x0=x0,x1=x0,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=cols[i])
      segments(x0=x1,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=cols[i])
    }
    if (missing(names.arg)) names.arg <- names(gene)
    if (!is.null(names.arg)) text(rep(xlim[1],length(y0)),y0+.5*(y0[2]-y0[1]),names(gene),cex=cex,adj=0)
  }
)

setMethod("genePlot",signature(gene='CompressedIRangesList'),
  function(gene, genome, refflat, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', ...) {
    if (missing(xlim)) xlim <- range(c(start(unlist(gene)),end(unlist(gene))))
    if (missing(names.arg)) names.arg <- names(gene)
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- length(gene)
    cols<-rep(rainbow(min(nsplice, 10)), ceiling(nsplice/10))
    y0 <- seq(1/(nsplice+1),1-1/(1+nsplice),by=1/(1+nsplice))
    for (i in 1:nsplice) {
      x0 <- start(gene[[i]])
      x1 <- end(gene[[i]])
      if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y=y0[i],y1=y0[i], col=cols[i])
      segments(x0=x0,x1=x1,y0=y0[i]+.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=cols[i])
      segments(x0=x0,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]-.33/(1+nsplice), col=cols[i])
      segments(x0=x0,x1=x0,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=cols[i])
      segments(x0=x1,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=cols[i])
    }
    if (missing(names.arg)) names.arg <- names(gene)
    if (!is.null(names.arg)) text(rep(xlim[1],length(y0)),y0+.5*(y0[2]-y0[1]),names.arg,cex=cex,adj=0)
  }
)


#getUCSCvariants2 <- function(gene, genome, refflat) {
#  if (missing(gene)) stop('gene must be specified')
#  if (length(gene)>1) stop('If argument gene is of type character, it must have length 1')
#  if (missing(refflat) & missing(genome)) stop('Either refflat or genome must be specified')
#  if (missing(refflat)) refflat <- getRefflat(genome=genome)
#  refflat <- refflat[refflat$geneName==gene,]
#  if (nrow(refflat)>0) {
#    exonStart <- strsplit(as.character(refflat$exonStarts),',')
#    exonEnd <- strsplit(as.character(refflat$exonEnds),',')
#    ans <- base::mapply(function(x,y) IRanges(as.numeric(x),as.numeric(y)), exonStart, exonEnd, SIMPLIFY=FALSE)
#    if (length(ans)) ans <- IRangesList(ans)
#  } else {
#    stop('Gene name not found in refflat table')
#  }
#  names(ans) <- as.character(refflat$name)
#  return(ans)
#}

setMethod("genePlot",signature(gene='character'),
  function(gene, genome, refflat, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', ...) {
    if (missing(gene)) stop('gene must be specified')
    if (length(gene)>1) stop('If argument gene is of type character, it must have length 1')
    if (missing(refflat) & missing(genome)) stop('Either refflat or genome must be specified')
    if (missing(refflat)) refflat <- getRefflat(genome=genome)
    refflat <- refflat[grep(gene,refflat$geneName),]
    if (nrow(refflat)>0) {
      exonStart <- strsplit(as.character(refflat$exonStarts),',')
      exonStart <- lapply(exonStart,function(x) as.numeric(x))
      exonEnd <- strsplit(as.character(refflat$exonEnds),',')
      exonEnd <- lapply(exonEnd,function(x) as.numeric(x))
      if (missing(xlim)) xlim <- range(c(exonStart,exonEnd))
      plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
      nsplice <- length(exonStart)
      y0 <- seq(1/(nsplice+1),1-1/(1+nsplice),by=1/(1+nsplice))
      for (i in 1:nsplice) {
        x0 <- exonStart[[i]]
        x1 <- exonEnd[[i]]
        segments(x0=x1[-length(x1)],x1=x0[-1],y=y0[i],y1=y0[i])
        segments(x0=x0,x1=x1,y0=y0[i]+.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice))
        segments(x0=x0,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]-.33/(1+nsplice))
        segments(x0=x0,x1=x0,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice))
        segments(x0=x1,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice))
      }
    } else {
      stop('Gene name not found in refflat table')
    }
  }
)


setMethod("rangesPlot",signature(x='IRangesList',gene='character'),
  function(x, gene, genome, refflat, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    if (length(gene)>1) stop('If argument gene is of type character, it must have length 1')
    if (missing(refflat) & missing(genome)) stop('Either refflat or genome must be specified')
    if (missing(refflat)) refflat <- getRefflat(genome=genome)
    refflat <- refflat[grep(gene,refflat$geneName),]
    if (nrow(refflat)>0) {
      txStart <- min(c(refflat$txStart,refflat$txEnd))
      txEnd <- max(c(refflat$txStart,refflat$txEnd))
      if (missing(xlim)) xlim <- c(txStart,txEnd)
      chr <- unique(as.character(refflat$chr))
      x <- x[[chr]]
      rangesPlot(x=x, gene=gene, refflat=refflat, exonProfile=exonProfile, xlab=xlab, ylab=ylab, xlim=xlim, heights=heights, ...)
    } else {
      stop('gene not found in specified genome/refflat')
    }
  }
)

setMethod("rangesPlot",signature(x='IRanges',gene='character'),
  function(x, gene, genome, refflat, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    if (length(gene)>1) stop('If argument gene is of type character, it must have length 1')
    if (missing(refflat) & missing(genome)) stop('Either refflat or genome must be specified')
    if (missing(refflat)) refflat <- getRefflat(genome=genome)
    refflat <- refflat[grep(gene,refflat$geneName),]
    if (nrow(refflat)>0) {
      txStart <- min(c(refflat$txStart,refflat$txEnd))
      txEnd <- max(c(refflat$txStart,refflat$txEnd))
      if (missing(xlim)) xlim <- c(txStart,txEnd)
      x <- x[start(x)>=xlim[1] & end(x)<=xlim[2]]
      x <- x[order(start(x))]
      cover <- coverage(x)
      cover <- seqselect(cover,(cover@lengths[1]+1):length(cover))
      layout(c(2,1),heights=heights)
      par(oma=c(1,0,0,0),mar=c(2.1,.1,.1,.1))
      genePlot(gene=gene, genome=genome, refflat=refflat, xlab='', ylab='', xlim=xlim, ...)
      if (exonProfile) {
        sel <- width(x)<=maxFragLength
        st <- c(start(x[sel]),start(x[!sel]),end(x[!sel])); en <- c(end(x[sel]),start(x[!sel]),end(x[!sel]))
        pr <- coverage(IRanges(start=st,end=en))
        pr <- seqselect(pr, (pr@lengths[1]+1):length(pr))
        minst <- min(st)
        lines(minst:(minst+length(pr)-1),pr/max(pr),col='gray')
      }
      #Height for large reads (i.e. spanning over several exons)
      x0 <- start(x); x1 <- end(x)
      y0 <- runif(sum(!sel))
      plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab='',ylab='',xaxt='n',yaxt='n',mgp=c(0,0,0))
      segments(x0=x0[!sel],x1=x1[!sel],y0=y0,y1=y0)
      #Height for short reads
      y0 <- as.integer(cover)[start(x)[sel]-min(start(x))+1] + runif(sum(sel))
      y0 <- y0/max(y0)
      segments(x0=x0[sel],x1=x1[sel],y0=y0,y1=y0,col='gray')
    } else {
      stop('gene not found in specified genome/refflat')
    }
  }
)


setMethod("rangesPlot",signature(x='IRanges',gene='IRanges'),
  function(x, gene, genome, refflat, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    txStart <- min(start(gene))
    txEnd <- max(end(gene))
    if (missing(xlim)) xlim <- c(txStart,txEnd)
    x <- x[start(x)>=xlim[1] & end(x)<=xlim[2]]
    if (length(x)==0) {
      warning('There are no reads in the specified region')
    } else {
      x <- x[order(start(x))]
      cover <- coverage(x)
      cover <- seqselect(cover,(cover@lengths[1]+1):length(cover))
      layout(c(2,1),heights=heights)
      par(oma=c(1,0,0,0),mar=c(2.1,.1,.1,.1))
      genePlot(gene=gene, xlab='', ylab='', xlim=xlim, ...)
      sel <- width(x)<=maxFragLength
      if (exonProfile) {
        st <- c(start(x[sel]),start(x[!sel]),end(x[!sel])); en <- c(end(x[sel]),start(x[!sel]),end(x[!sel]))
        pr <- coverage(IRanges(start=st,end=en))
        pr <- seqselect(pr, (pr@lengths[1]+1):length(pr))
        minst <- min(st)
        lines(minst:(minst+length(pr)-1),pr/max(pr),col='gray')
      }
      #Height for large reads (i.e. spanning over several exons)
      plot(NA,NA,xlim=xlim,ylim=0:1,xlab='',ylab='',xaxt='n',yaxt='n',mgp=c(0,0,0))
      x0 <- start(x); x1 <- end(x)
      if (any(!sel)) {
        y0 <- runif(sum(!sel))
        segments(x0=x0[!sel],x1=x1[!sel],y0=y0,y1=y0)
      }
      #Height for short reads
      if (any(sel)) {
        y0 <- as.integer(cover)[start(x)[sel]-start(x)[1]+1] + runif(sum(sel))
        y0 <- y0/max(y0)
        segments(x0=x0[sel],x1=x1[sel],y0=y0,y1=y0,col='gray')
      }
    }
  }
)


setMethod("rangesPlot",signature(x='IRanges',gene='IRangesList'),
  function(x, gene, genome, refflat, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    txStart <- min(unlist(start(gene)))
    txEnd <- max(unlist(end(gene)))
    if (missing(xlim)) xlim <- c(txStart,txEnd)
    x <- x[start(x)>=xlim[1] & end(x)<=xlim[2]]
    if (length(x)==0) {
      warning('There are no reads in the specified region')
    } else {
      x <- x[order(start(x))]
      cover <- coverage(x)
      cover <- seqselect(cover,(cover@lengths[1]+1):length(cover))
      layout(c(2,1),heights=heights)
      par(oma=c(1,0,0,0),mar=c(2.1,.1,.1,.1))
      genePlot(gene=gene, xlab='', ylab='', xlim=xlim, ...)
      sel <- width(x)<=maxFragLength
      if (exonProfile) {
        st <- c(start(x[sel]),start(x[!sel]),end(x[!sel])); en <- c(end(x[sel]),start(x[!sel]),end(x[!sel]))
        pr <- coverage(IRanges(start=st,end=en))
        pr <- seqselect(pr, (pr@lengths[1]+1):length(pr))
        minst <- min(st)
        lines(minst:(minst+length(pr)-1),pr/max(pr),col='gray')
      }
      #Height for large reads (i.e. spanning over several exons)
      plot(NA,NA,xlim=xlim,ylim=0:1,xlab='',ylab='',xaxt='n',yaxt='n',mgp=c(0,0,0))
      x0 <- start(x); x1 <- end(x)
      if (any(!sel)) {
        y0 <- as.integer(cover)[start(x)[!sel]-min(start(x))+1] + runif(sum(!sel))
        y0 <- y0/max(y0)
       # y0 <- runif(sum(!sel))
        segments(x0=x0[!sel],x1=x1[!sel],y0=y0,y1=y0)
      }
      #Height for short reads
      if (any(sel)) {
        y0 <- as.integer(cover)[start(x)[sel]-min(start(x))+1] + runif(sum(sel))
        y0 <- y0/max(y0)
        segments(x0=x0[sel],x1=x1[sel],y0=y0,y1=y0,col='gray')
      }
    }
  }
)


setMethod("rangesPlot",signature(x='RangedData',gene='character'),
  function(x, gene, genome, refflat, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    rangesPlot(unlist(ranges(x)),gene=gene,genome=genome,refflat=refflat,exonProfile=exonProfile,maxFragLength=maxFragLength,xlab=xlab,ylab=ylab,xlim=xlim,heights=heights,...)
  }
)

setMethod("rangesPlot",signature(x='RangedData',gene='IRanges'),
  function(x, gene, genome, refflat, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    rangesPlot(unlist(ranges(x)),gene=gene,genome=genome,refflat=refflat,exonProfile=exonProfile,maxFragLength=maxFragLength,xlab=xlab,ylab=ylab,xlim=xlim,heights=heights,...) 
  }
)

setMethod("rangesPlot",signature(x='RangedData',gene='IRangesList'),
function(x, gene, genome, refflat, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
  rangesPlot(unlist(ranges(x)),gene=gene,genome=genome,refflat=refflat,exonProfile=exonProfile,maxFragLength=maxFragLength,xlab=xlab,ylab=ylab,xlim=xlim,heights=heights,...)  
}
)
