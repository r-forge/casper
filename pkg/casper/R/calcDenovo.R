require(methods)

#Class denovoGeneExpr
setClass("denovoGeneExpr", representation(posprob = "matrix", expression = "data.frame", variants = "RangedData"))

valid_denovoGeneExpr <- function(object) {
  msg <- NULL
  if (any(!(c('model','posprob') %in% colnames(object@posprob)))) msg <- "Incorrect colnames in 'posprob'"
  if (any(!(c('model','expr','varName') %in% colnames(object@expression)))) msg <- "Incorrect colnames in 'expression'"
  if (class(object@variants)!='RangedData') msg <- "Element variants must be of class 'RangedData'"
  if (is.null(msg)) { TRUE } else { msg }
}

setValidity("denovoGeneExpr", valid_denovoGeneExpr)

setMethod("show", signature(object="denovoGeneExpr"), function(object) {
  cat("denovoGeneExpr object\n")
  cat("\nPosterior model probabilities\n")
  show(object@posprob)
  cat("\nEstimated expression (conditional on each model)\n")
  show(object@expression)
  cat("\nUse 'variants' method to access exon starts/ends for all variants\n")
}
)

setGeneric("variants", function(object) standardGeneric("variants"))
setMethod("variants", signature(object="denovoGeneExpr"), function(object) {
  object@variants
}
)

          
#Class denovoGenomeExpr
setClass("denovoGenomeExpr", representation(genes = "list"))

valid_denovoGenomeExpr <- function(object) {
  msg <- NULL
  if (any(sapply(object@genes,class)!='denovoGeneExpr')) msg <- "All elements must be of class denovoGeneExpr"
  if (is.null(msg)) { TRUE } else { msg }
}

setValidity("denovoGenomeExpr", valid_denovoGenomeExpr)

setMethod("show", signature(object="denovoGenomeExpr"), function(object) {
  cat("denovoGenomeExpr object with",length(object@genes),"gene islands\n")
}
)


calcDenovo <- function(distrs, genomeDB, pc, readLength, geneid, priorq=3, mprior, minpp=0.01, selectBest=FALSE, verbose=FALSE, mc.cores=1) {
  if (missing(readLength)) stop("readLength must be specified")
  if (class(genomeDB)!='denovoGenome') stop("genomeDB must be of class 'denovoGenome'")
  if (!all(c('nvarPrior','nexonPrior') %in% names(mprior))) stop("Incorrect mprior. Please use modelPrior to generate it.")
  
  #Format input
  startcdf <- as.double(ecdf(distrs$stDis)(seq(0,1,.001)))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  nvarPrior <- mprior$nvarPrior$nbpar  #format as list 1:40, larger numbers just use 40th elem. Ensure all 1:40 are included.
  nexonPrior <- mprior$nexonPrior$bbpar
  minpp <- as.double(minpp)
  selectBest <- as.integer(selectBest)
  verbose <- as.integer(verbose)
  if (missing(geneid)) geneid <- names(genomeDB@genes)[sapply(genomeDB@genes,length)>1]

  exons <- lapply(genomeDB@genes,function(z) as.integer(names(z)))
  exonwidth <- lapply(genomeDB@genes,width)

  if (!all(geneid %in% names(exons))) stop('geneid not found in genomeDB@exons')
  if (!all(geneid %in% names(pc))) stop('geneid not found in pc')
  if (!all(geneid %in% names(genomeDB@transcripts))) stop('geneid not found in genomeDB@transcripts')

  #Define basic function
  f <- function(z) {
    geneid <- as.integer(z)
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- genomeDB@transcripts[z]
    pc <- pc[z]
    ans <- calcDenovoMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,geneid=as.list(geneid),pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,nvarPrior=nvarPrior,nexonPrior=nexonPrior,priorq=priorq,minpp=minpp,selectBest=selectBest,verbose=verbose)
    mapply(function(z1,z2) formatDenovoOut(z1,z2), ans, genomeDB@genes[z], SIMPLIFY=FALSE)
  }

  #Run
  if (mc.cores>1) {
    if ('multicore' %in% loadedNamespaces()) {
      #split into smaller jobs
      nsplit <- floor(length(geneid)/mc.cores)
      geneid <- lapply(1:mc.cores, function(z) { geneid[((z-1)*nsplit+1):min((z*nsplit),length(geneid))] })
      ans <- mclapply(geneid,f,mc.cores=mc.cores)
      ans <- do.call(c,ans)
     } else stop('multicore library has not been loaded!')
  } else {
    ans <- f(geneid)
    names(ans) <- geneid
  }
  new("denovoGenomeExpr", genes=ans)
}



formatDenovoOut <- function(ans, genesel) {
  colnames(ans[[1]]) <- c('model','posprob')
  ans[[2]] <- data.frame(ans[[2]],ans[[3]])
  ans[[3]] <- NULL
  colnames(ans[[2]]) <- c('model','expr','varName')
  variants <- genesel[as.character(ans[[3]])]
  names(variants) <- NULL
  ans[[3]] <- RangedData(ranges=variants,space=ans[[4]])
  ans[[4]] <- NULL
  names(ans) <- c('posprob','expression','variants')
  new("denovoGeneExpr",posprob=ans$posprob,expression=ans$expression,variants=ans$variants)
}

#calcDenovoSingle <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, nvarPrior, nexonPrior, priorq, minpp, selectBest, verbose) {
#  ans <- .Call("calcDenovoSingle",as.integer(exons),as.integer(exonwidth),transcripts,geneid,pc,startcdf,lendis,lenvals,readLength,nvarPrior,nexonPrior,priorq,minpp,selectBest,verbose)
#  return(ans)
#}

calcDenovoMultiple <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, nvarPrior, nexonPrior, priorq, minpp, selectBest, verbose) {
  ans <- .Call("calcDenovoMultiple",exons,exonwidth,transcripts,geneid,pc,startcdf,lendis,lenvals,readLength,nvarPrior,nexonPrior,priorq,minpp,selectBest,verbose)
  return(ans)
}


