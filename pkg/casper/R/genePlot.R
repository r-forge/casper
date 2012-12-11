buildGene<-function(genomeDB, islandid, txs){
    rd<-lapply(txs, function(x) genomeDB@islands[[islandid]][names(genomeDB@islands[[islandid]]) %in% as.character(x),])
    rdl<-IRangesList(rd)
    rdl
}

setMethod("genePlot", signature(generanges='missing',islandid='character',genomeDB='annotatedGenome',reads='missing',exp='missing') , function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
  txs <- transcripts(islandid=islandid,genomeDB=genomeDB)
  genePlot(txs, names.arg=names.arg, xlab=xlab, ylab=ylab, xlim=xlim, cex=cex, yaxt=yaxt, col=col, ...)
}
)

setMethod("genePlot", signature(generanges='missing',islandid='character',genomeDB='annotatedGenome',reads='procBam',exp='ExpressionSet') , function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
  txexp<-exprs(exp)[names(genomeDB@transcripts[[islandid]]),]
  if(sum(txexp)>0) txexp<-(txexp/sum(txexp))*100
  gene<-buildGene(txs=genomeDB@transcripts[[islandid]], genomeDB, islandid=islandid)
  chr <- genomeDB@exon2island[match(TRUE,genomeDB@exon2island==islandid),'space']
  reads@pbam <- reads@pbam[chr]
  rangesPlot(x=reads, gene, exonProfile=FALSE)
}
)

setMethod("genePlot", signature(generanges='missing',islandid='character',genomeDB='annotatedGenome',reads='RangedData',exp='ExpressionSet') , function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
  txexp<-exprs(exp)[names(genomeDB@transcripts[[islandid]]),]
  if(sum(txexp)>0) txexp<-(txexp/sum(txexp))*100
  gene<-buildGene(txs=genomeDB@transcripts[[islandid]], genomeDB, islandid=islandid)
  rangesPlot(x=reads[as.character(genomeDB@exon2island[genomeDB@exon2island$id==names(gene[[1]])[1],]$space)], gene, exonProfile=FALSE)
  #txexp
}
)

setMethod("plotExpr",signature(gene="denovoGeneExpr"),
  function(gene, minProbExpr, minExpr, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    n <- names(variants(gene))
    n[nchar(n)>30] <- paste("Variant ",1:sum(nchar(n)>30),sep="")
    names(n) <- names(variants(gene))
    v <- variantMargExpr(gene,minProbExpr=minProbExpr,minExpr=minExpr)
    names.arg <- n[rownames(v)]
    names.arg <- paste(names.arg,' (Rel. expr.=',round(v[,'expr'],3),'; P(expressed)=',round(v[,'probExpressed'],3),')',sep='')
    names(names.arg) <- n[rownames(v)]
    gene <- variants(gene)[rownames(v)]
    o <- names(n)[names(n) %in% names(gene)]
    gene <- gene[o]
    names.arg <- names.arg[n[o]]
    start(gene) <- round(start(gene)/1000)
    end(gene) <- round(end(gene)/1000)
    if (missing(xlim)) xlim <- c(min(start(gene)),max(end(gene)))
    genePlot(gene, names.arg=names.arg, xlab=xlab, ylab=ylab, xlim=xlim, cex=1, yaxt='n', col, ...)
  }
)


setMethod("genePlot",signature(generanges='IRanges'),
  function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    if (missing(xlim)) xlim <- range(c(start(generanges),end(generanges)))
    if (missing(col)) col <- 1
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- 1
    y0 <- seq(1/(1+nsplice),1-1/(1+nsplice),by=1/(1+nsplice))
    x0 <- start(generanges)
    x1 <- end(generanges)
    if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y0=y0,y1=y0)
    segments(x0=x0,x1=x1,y0=y0+.33/(1+nsplice),y1=y0+.33/(1+nsplice),col=col)
    segments(x0=x0,x1=x1,y0=y0-.33/(1+nsplice),y1=y0-.33/(1+nsplice),col=col)
    segments(x0=x0,x1=x0,y0=y0-.33/(1+nsplice),y1=y0+.33/(1+nsplice),col=col)
    segments(x0=x1,x1=x1,y0=y0-.33/(1+nsplice),y1=y0+.33/(1+nsplice),col=col)
  }
)


setMethod("genePlot",signature(generanges='RangedData'),
  function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    genePlot(ranges(generanges), names.arg=names.arg, xlab=xlab, ylab=ylab, xlim=xlim, cex=cex, yaxt=yaxt, col=col, ...)
  }
)

setMethod("genePlot",signature(generanges='IRangesList'),
  function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    if (missing(xlim)) xlim <- range(c(start(generanges),end(generanges)))
    if (missing(names.arg)) names.arg <- names(generanges)
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- length(generanges)
    if (missing(col)) col<-rep(rainbow(min(nsplice, 10)), ceiling(nsplice/10)) else if (length(col)!=nsplice) col <- rep(col[1],nsplice)
    y0 <- seq(1/(1+nsplice),1-1/(1+nsplice),by=1/(1+nsplice))
    for (i in nsplice:1) {
      x0 <- start(generanges[[i]])
      x1 <- end(generanges[[i]])
      if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y0=y0[i],y1=y0[i], col=col[i])
      segments(x0=x0,x1=x1,y0=y0[i]+.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
      segments(x0=x0,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]-.33/(1+nsplice), col=col[i])
      segments(x0=x0,x1=x0,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
      segments(x0=x1,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
    }
    if (!is.null(names.arg)) text(rep(xlim[1],length(y0)),y0+.5*(y0[2]-y0[1]),names.arg,cex=cex,adj=0)
  }
)

setMethod("genePlot",signature(generanges='CompressedIRangesList'),
  function(generanges, islandid, genomeDB, reads, exp, names.arg, xlab='', ylab='', xlim, cex=1, yaxt='n', col, ...) {
    if (missing(xlim)) xlim <- range(c(start(unlist(generanges)),end(unlist(generanges))))
    if (missing(names.arg)) names.arg <- names(generanges)
    plot(NA,NA,xlim=xlim,ylim=c(0,1),xlab=xlab,ylab=ylab,yaxt=yaxt,...)
    nsplice <- length(generanges)
    if (missing(col)) col<-rep(rainbow(min(nsplice, 10)), ceiling(nsplice/10)) else if (length(col)!=nsplice) col <- rep(col[1],nsplice)
    y0 <- rev(seq(1/(1+nsplice),1-1/(1+nsplice),by=1/(1+nsplice)))
    for (i in 1:nsplice) {
      x0 <- start(generanges[[i]])
      x1 <- end(generanges[[i]])
      if (length(x0)>1) segments(x0=x1[-length(x1)],x1=x0[-1],y0=y0[i],y1=y0[i], col=col[i])
      segments(x0=x0,x1=x1,y0=y0[i]+.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
      segments(x0=x0,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]-.33/(1+nsplice), col=col[i])
      segments(x0=x0,x1=x0,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
      segments(x0=x1,x1=x1,y0=y0[i]-.33/(1+nsplice),y1=y0[i]+.33/(1+nsplice), col=col[i])
    }
    if (!is.null(names.arg)) text(rep(xlim[1],length(y0)),y0+.5*(y0[1]-y0[2]),names.arg,cex=cex,adj=0)
  }
)




#setMethod("rangesPlot",signature(x='IRanges',gene='IRanges'),
#  function(x, gene, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
#    txStart <- min(start(gene))
#    txEnd <- max(end(gene))
#    if (missing(xlim)) xlim <- c(txStart,txEnd)
#    x <- x[start(x)>=xlim[1] & end(x)<=xlim[2]]
#    if (length(x)==0) {
#      warning('There are no reads in the specified region')
#    } else {
#      x <- x[order(start(x))]
#      cover <- coverage(x)
#      cover <- seqselect(cover,(cover@lengths[1]+1):length(cover))
#      layout(c(2,1),heights=heights)
#      par(oma=c(1,0,0,0),mar=c(2.1,.1,.1,.1))
#      genePlot(gene=gene, xlab='', ylab='', xlim=xlim, ...)
#      sel <- width(x)<=maxFragLength
#      if (exonProfile) {
#        st <- c(start(x[sel]),start(x[!sel]),end(x[!sel])); en <- c(end(x[sel]),start(x[!sel]),end(x[!sel]))
#        pr <- coverage(IRanges(start=st,end=en))
#        pr <- seqselect(pr, (pr@lengths[1]+1):length(pr))
#        minst <- min(st)
#        lines(minst:(minst+length(pr)-1),pr/max(pr),col='gray')
#      }
#      #Height for large reads (i.e. spanning over several exons)
#      plot(NA,NA,xlim=xlim,ylim=0:1,xlab='',ylab='',xaxt='n',yaxt='n',mgp=c(0,0,0))
#      x0 <- start(x); x1 <- end(x)
#      if (any(!sel)) {
#        y0 <- runif(sum(!sel))
#        segments(x0=x0[!sel],x1=x1[!sel],y0=y0,y1=y0)
#      }
#      #Height for short reads
#      if (any(sel)) {
#        y0 <- as.integer(cover)[start(x)[sel]-start(x)[1]+1] + runif(sum(sel))
#        y0 <- y0/max(y0)
#        segments(x0=x0[sel],x1=x1[sel],y0=y0,y1=y0,col='gray')
#      }
#    }
#  }
#)


setMethod("rangesPlot",signature(x='IRanges'),
  function(x, gene, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
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


setMethod("rangesPlot",signature(x='RangedData',gene='IRanges'),
  function(x, gene, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    if (missing(xlim)) xlim <- c(min(start(gene)),max(end(gene)))
    rangesPlot(unlist(ranges(x)),gene=gene,exonProfile=exonProfile,maxFragLength=maxFragLength,xlab=xlab,ylab=ylab,xlim=xlim,heights=heights,...) 
  }
)

setMethod("rangesPlot",signature(x='RangedData',gene='IRangesList'),
function(x, gene, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
  if (missing(xlim)) xlim <- c(min(unlist(start(gene))),max(unlist(end(gene))))
  rangesPlot(unlist(ranges(x)),gene=gene,exonProfile=exonProfile,maxFragLength=maxFragLength,xlab=xlab,ylab=ylab,xlim=xlim,heights=heights,...)  
}
)



setMethod("rangesPlot",signature(x='procBam'),
  function(x, gene, exonProfile=TRUE, maxFragLength=300, xlab='', ylab='', xlim, heights=c(2,1), ...) {
    x <- getReads(x)
    x <- x[start(x)>=xlim[1] & end(x)<=xlim[2],]
    if (nrow(x)==0) {
      warning('There are no reads in the specified region')
    } else {
      layout(c(2,1),heights=heights)
      par(oma=c(1,0,0,0),mar=c(2.1,.1,.1,.1))
      genePlot(gene=gene, xlab='', ylab='', xlim=xlim, ...)
      tab <- table(x[['id']])
      sel <- x[['id']] %in% c(names(tab[tab>2]), findLongInserts(x, minsize=500))
      notsel <- !sel
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
      y0 <- x[['id']]/max(x[['id']])
      if (any(notsel)) {
        segments(x0=x0[notsel],x1=x1[notsel],y0=y0[notsel])
      }
      #Height for short reads
      if (any(sel)) {
        segments(x0=x0[sel],x1=x1[sel],y0=y0[sel],col='gray')
      }
    }
  }
)



findLongInserts <- function(x, minsize) {
  #Select long paired end inserts in a RangedData object ordered by pair id (e.g. as returned by procBam in package casper)
  # - x: RangedData, column id must indicate read pair id. Assumed to be ordered according to read id & fragment start
  # - minsize: minimum length for an insert to be considered long. Defaults to 99% percentile of observed insert widths
  # Returns: ids of long inserts
  w <- end(x)[-1] - start(x)[-nrow(x)]
  sel <- x$id[-nrow(x)] == x$id[-1]
  w[!sel] <- 0
  if (missing(minsize)) minsize <- quantile(w[sel],probs=.99)
  #Select long inserts
  idsel <- x$id[which(w>minsize)]
  sel <- which(x$id %in% idsel)
  chr <- x$space[sel]; st <- start(x)[sel]; en <- end(x)[sel]; id <- x$id[sel]
  leftend <- which(c(TRUE,id[-1]!=id[-length(id)]))
  rightend <- c(leftend[-1]-1,length(id))
  st <- st[leftend]; en <- en[rightend]; chr <- chr[leftend]; id <- id[leftend]
  id[st<en]
}
