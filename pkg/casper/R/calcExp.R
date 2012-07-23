proc <- function(distrs, genomeDB, pc, readLength, geneid, relativeExpr=TRUE, priorq=3, report=0, niter=10^3, burnin=100, mc.cores=1) {
  #Format input
  startcdf <- distrs$stDis(seq(0,1,.001))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  report <- as.integer(report)
  niter <- as.integer(niter)
  burnin <- as.integer(burnin)
  if (!report %in% 0:2) stop("Argument report must be equal to 0, 1 or 2")
  if ((niter<=burnin) & report>=2) stop("Too many burnin iterations specified")
  if (missing(geneid)) geneid <- names(genomeDB@islands)[sapply(genomeDB@islands,length)>1]
  
  exons <- lapply(genomeDB@islands,function(z) as.integer(names(z)))
  exonwidth <- lapply(genomeDB@islands,width)
  strand <- genomeDB@islandStrand

  if (!all(geneid %in% names(exons))) stop('geneid not found in genomeDB@islands')
  if (!all(geneid %in% names(pc))) stop('geneid not found in pc')
  if (!all(geneid %in% names(genomeDB@transcripts))) stop('geneid not found in genomeDB@transcripts')
  if (!all(geneid %in% names(genomeDB@islandStrand))) stop('geneid not found in genomeDB@islandStrand')
  

      #Format input
      startcdf <- distrs$stDis(seq(0,1,.001))
      lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
      lenvals <- as.integer(names(distrs$lenDis))
      readLength <- as.integer(readLength)
      priorq <- as.double(priorq)
      citype <- as.integer(switch(citype, none=0, asymp=1, exact=2))
      if (length(citype)==0) stop("citype should take the value 'none', 'asymp', or 'exact'")
      niter <- as.integer(niter)
      burnin <- as.integer(burnin)
      verbose <- as.integer(verbose)
      if ((niter<=burnin) & citype==2) stop("Too many burnin iterations specified. Decrease burnin or increase niter")
      if (missing(geneid)) geneid <- names(genomeDB@islands)[sapply(genomeDB@islands,length)>1]

      exons <- lapply(genomeDB@islands,function(z) as.integer(names(z)))
      exonwidth <- lapply(genomeDB@islands,width)
      strand <- genomeDB@islandStrand

      if (!all(geneid %in% names(exons))) stop('geneid not found in genomeDB@islands')
      if (!all(geneid %in% names(pc@counts))) stop('geneid not found in pc')
      if (!all(geneid %in% names(genomeDB@transcripts))) stop('geneid not found in genomeDB@transcripts')
      if (!all(geneid %in% names(genomeDB@islandStrand))) stop('geneid not found in genomeDB@islandStrand')


      #Define basic function

  f <- function(z) {
    geneid <- as.integer(z)
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- genomeDB@transcripts[z]
    strand <- as.list(as.integer(ifelse(strand[z]=='+', 1, -1)))
    pc <- pc[z]
    ans <- calcKnownMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,geneid=as.list(geneid),pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,priorq=priorq, strand=strand, report=report, niter=niter, burnin=burnin)
    if (report==0) {
      ans <- lapply(ans, function(z) { res=vector("list",1); res[[1]]= z[[1]]; names(res[[1]])= z[[2]]; res })
    } else if (report==1) {
      ans <- lapply(ans, function(z) { res=vector("list",2); res[[1]]= z[[1]]; res[[2]]= z[[3]]; names(res[[1]])= names(res[[2]])= z[[2]]; res })
    } else if (report==2) {
      ans <- lapply(ans, function(z) { res=vector("list",3); res[[1]]= z[[1]]; res[[2]]= z[[3]]; res[[3]]= matrix(z[[4]],nrow=niter-burnin); names(res[[1]])= names(res[[2]])= z[[2]]; res })
    }
    ans
  }


    #Run

  sel <- !sapply(pc[geneid], is.null)
  all <- geneid
  geneid <- geneid[sel]

    sel <- !sapply(pc@counts[geneid], is.null)
    all <- geneid
    geneid <- geneid[sel]
    if (verbose) cat("Obtaining expression estimates...\n")
    if (mc.cores>1 && length(geneid)>mc.cores) {
      if ('multicore' %in% loadedNamespaces()) {
                    #split into smaller jobs
        nsplit <- floor(length(geneid)/mc.cores)
        geneid <- lapply(1:mc.cores, function(z) { geneid[((z-1)*nsplit+1):min((z*nsplit),length(geneid))] })
        ans <- mclapply(geneid,f,mc.cores=mc.cores)
        ans <- unlist(ans, recursive=F)
        names(ans) <- unlist(geneid)
      } else stop('multicore library has not been loaded!')
    } else {
      ans <- f(geneid)
      names(ans) <- geneid
    }
    miss <- lapply(genomeDB@transcripts[all[!(all %in% names(ans))]], names)

  
    #Format as ExpressionSet

    if (verbose) cat("Formatting output...\n")
    transcript <- c(unlist(lapply(ans, function(z) names(z[[1]]))), unname(unlist(miss)))
    gene <- c(rep(names(ans),sapply(ans,function(z) length(z[[1]]))), rep(names(miss), sapply(miss, length)))
    fdata <- data.frame(transcript=transcript, gene=gene)
    misse <- rep(NA, length(unlist(miss)))
    names(misse) <- unlist(miss)
    exprsx <- matrix(c(unlist(lapply(ans,'[[',1)), misse) ,ncol=1)
    rownames(exprsx) <- c(unlist(lapply(ans, function(z) names(z[[1]]))), names(misse))
    exprsx[exprsx==-1] <- NA
    if (citype>=1) {
      se <- c(unlist(lapply(ans,'[[',2)), misse)
      se[se==-1] <- NA
      fdata$ci95.low <- exprsx - 1.96*se; fdata$ci95.low[fdata$ci95.low<0] <- 0
      fdata$ci95.high <- exprsx + 1.96*se; fdata$ci95.high[fdata$ci95.high>1] <- 1
    }
    if (citype>=2) {
      q <- lapply(ans,function(z) apply(z[[3]],2,quantile,probs=c(.025,.975)))
      q[names(miss)] <- c(NA, NA)
      q <- t(do.call(cbind,q))
      fdata$ci95.low <- q[,1]; fdata$ci95.high <- q[,2]
    }
    rownames(exprsx) <- rownames(fdata) <- fdata$transcript
    fdata <- new("AnnotatedDataFrame",fdata)
    ans <- new("ExpressionSet",exprs=exprsx,featureData=fdata)

    #Return absolute expression levels
  if (!relativeExpr) {
    nreads <- sapply(pc[unique(gene)],sum)
    exprs(ans) <- exprs(ans)*nreads[ans@featureData$gene]
    exprs(ans) <- t(t(exprs(ans))/colSums(exprs(ans), na.rm=T))
  }
  ans
}

genomeBystrand <- function(DB, strand){
  sel <- DB@islandStrand==strand
  islands <- DB@islands[sel]
  transcripts <- DB@transcripts[sel]
  exonsNI <- DB@exonsNI[DB@exonsNI$id %in% unlist(transcripts),]
  exon2island <- DB@exon2island[DB@exon2island$id %in%unlist(transcripts),]
  islandStrand <- DB@islandStrand[sel]
  txid <- unlist(lapply(transcripts, names))
  aliases <- DB@aliases[DB@aliases$tx %in% txid,]
  ans <- new("annotatedGenome", aliases=aliases, denovo=TRUE, exonsNI=exonsNI, islandStrand=islandStrand, transcripts=transcripts, exon2island=exon2island, dateCreated=Sys.Date(), genomeVersion=DB@genomeVersion, islands=islands)
  ans
}

mergeExp <- function(minus, plus){
  exp <- rbind(exprs(plus), exprs(minus))
  fdata <- new("AnnotatedDataFrame", rbind(fData(plus), fData(minus)))
  ans <- new("ExpressionSet",exprs=exp,featureData=fdata)
  ans
}

calcExp <- function(distrs, genomeDB, pc, readLength, geneid, relativeExpr=TRUE, priorq=3, report=0, niter=10^3, burnin=100, mc.cores=1) {
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if (pc@denovo) stop("pc must be a pathCounts object from known genome")
  if(pc@stranded){
    plusDB <- genomeBystrand(genomeDB, "+")
    plusGI <- geneid[geneid %in% names(plusDB@transcripts)]
    plus <- proc(distrs, plusDB, pc@counts$plus, readLength=readLength, geneid=plusGI, relativeExpr=relativeExpr, priorq=priorq, report=report, niter=niter, burnin=burnin, mc.cores=mc.cores)
    minusDB <- genomeBystrand(genomeDB, "-")
    minusGI <- geneid[geneid %in% names(minusDB@transcripts)]
    minus <- proc(distrs, minusDB, pc=pc@counts$minus, readLength=readLength, minusGI, relativeExpr=relativeExpr, priorq=priorq, report=report, niter=niter, burnin=burnin, mc.cores=mc.cores)
    ans <- mergeExp(plus, minus)
  } else {
    ans <- proc(distrs, genomeDB, pc@counts[[1]], readLength=readLength, geneid, relativeExpr=relativeExpr, priorq=priorq, report=report, niter=niter, burnin=burnin, mc.cores=mc.cores)
  }
  ans
}
    
calcKnownMultiple <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, priorq, strand, citype, niter, burnin, verbose) {
  ans <- .Call("calcKnownMultiple",exons,exonwidth,transcripts,geneid,pc,startcdf, lendis, lenvals, readLength, priorq, strand, citype, niter, burnin, verbose)
  return(ans)
}


lhoodGrid <- function(pc, distrs, genomeDB, readLength, geneid, grid, priorq=2) {
  if (length(geneid)>1) stop("Only 1 gene is allowed, please specify a single gene in argument 'geneid'")
  if (missing(readLength)) stop("readLength must be specified")
  if (genomeDB@denovo) stop("genomeDB must be a known genome")
  if (pc@denovo) stop("pc must be a pathCounts object from known genome")

  #Format input
  startcdf <- distrs$stDis(seq(0,1,.001))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)

  exons <- as.integer(names(genomeDB@islands[[geneid]]))
  exonwidth <- width(genomeDB@islands[[geneid]])
  strand <- genomeDB@islandStrand[geneid]
  geneid <- as.integer(geneid)
  transcripts <- genomeDB@transcripts[[geneid]]
  if (length(transcripts)==1) stop("Single transcript specified, estimation not run")
  strand <- as.integer(ifelse(strand=='+', 1, -1))
  pc <- pc@counts[[geneid]]
  
  if (missing(grid)) {
    grid <- lapply(1:(length(transcripts)-1), function(z) seq(-10,10,length=10^(5/(length(transcripts)-1))))
    #grid <- lapply(1:length(transcripts), function(z) seq(0.0001,.9999,length=10^(5/length(transcripts))))
  } else {
    if (!is.list(grid)) stop("grid must be of class list")
    if (length(grid) != length(transcripts)-1) stop("grid must have one row less than the number of known transcripts")
  }
  milogit <- function(theta) cbind(rep(1,nrow(theta)),exp(theta))/(1+rowSums(exp(theta))) #matrix conversion
  gridmat <- t(milogit(expand.grid(grid)))
  
  ans <- .Call("lhoodGrid",gridmat,exons,exonwidth,transcripts,pc,startcdf,lendis,lenvals,readLength,priorq,strand)
  names(ans[[2]]) <- ans[[3]]; ans[[3]] <- NULL
  names(ans) <- c('logpos','piest','logpos.piest','S','pathprob','marginalLhood')
  rownames(ans$pathprob) <- names(ans$piest)
  ans$pathprob <- t(ans$pathprob)
  names(ans$marginalLhood) <- c('laplace','IS')
  ans$grid <- grid
  return(ans)  
}


calcExpOld<-function(distrs, genomeDB, pc, readLength){
    if (missing(readLength)) stop("readLength must be specified")
      if (class(genomeDB)!="knownGenome") stop("genomeDB must be of class 'knownGenome'")
      dexo<-as.data.frame(genomeDB@exonsNI)
      mexo<-as.matrix(dexo[,c(5,4)])
      startcdf <- as.double(ecdf(distrs$stDis)(seq(0,1,.001)))
      #stafun = ecdf(distrs$stDis)
      #fill<-c(rep(0,(min(as.numeric(names(distrs$lenDis))))), distrs$lenDis)
      #names(fill)[1:min(as.numeric(names(distrs$lenDis)))]<-0:(min(as.numeric(names(distrs$lenDis)))-1)
      lendis<- as.double(distrs$lenDis/sum(distrs$lenDis))
      exp<-.Call("calc", mexo, genomeDB@newTx, pc, startcdf, lendis, as.integer(names(lendis)), as.integer(readLength))
      names(exp)<-names(genomeDB@newTx)
      exp<-as.matrix(exp)
      featureData<-as.numeric(unlist(genomeDB@txs@elementMetadata$gene <- id))
      names(featureData)<-unlist(genomeDB@txs@elementMetadata$tx <- name)
      featureData<-as.data.frame(as.matrix(featureData[rownames(exp)]))
      metadata <- data.frame(labelDescription = "gene_id", row.names="gene_id")
      featureData<-new("AnnotatedDataFrame", data=featureData, varMetadata=metadata)
      exp<-new("ExpressionSet", exprs=exp, featureData=featureData, annotation="refseq")
      names(fData(exp)) <- 'entrezid'
      #Compute rpkm
      exons <- unlist(genomeDB@exons)
      geneid <- cumsum(values(exons)$exon <- rank==1)
      txlength <- tapply(width(exons),geneid,FUN=sum)
      names(txlength) <- names(genomeDB@exons)
      fData(exp)$txlength <- txlength[featureNames(exp)]
      rpkm <- matrix(10^9 * exprs(exp)[,] / fData(exp)$txlength, nrow=nrow(exp), ncol=ncol(exp))
      exp<-new("ExpressionSet", exprs=rpkm, featureData=featureData(exp), annotation="refseq")
      exp
  }

