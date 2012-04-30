#########################################################################
## DEFINITION AND METHODS FOR CLASS denovoGeneExpr AND denovoGenomeExpr
#########################################################################

require(methods)

setClass("denovoGeneExpr", representation(posprob = "data.frame", expression = "data.frame", variants = "RangedData", integralSum= "numeric", npathDeleted= "numeric"))

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
setClass("denovoGenomeExpr", representation(islands = "list"))

valid_denovoGenomeExpr <- function(object) {
  msg <- NULL
  if (any(sapply(object@islands,class)!='denovoGeneExpr')) msg <- "All elements must be of class denovoGeneExpr"
  if (is.null(msg)) { TRUE } else { msg }
}

setValidity("denovoGenomeExpr", valid_denovoGenomeExpr)

setMethod("show", signature(object="denovoGenomeExpr"), function(object) {
  cat("denovoGenomeExpr object with",length(object@islands),"gene islands\n")
}
)

setMethod("[", signature(x="denovoGenomeExpr"), function(x, i, ...) { new("denovoGenomeExpr", islands=x@islands[i]) })
setMethod("[[", signature(x="denovoGenomeExpr"), function(x, i, j, ...) { x@islands[[i]] } )

                            

#########################################################################
## Function calcDenovo
#########################################################################

calcDenovo <- function(distrs, genomeDB, pc, readLength, geneid, priorq=3, mprior, minpp=0.01, selectBest=FALSE, method='auto', niter=10^4, exactMarginal=TRUE, verbose=FALSE, mc.cores=1) {
  if (missing(readLength)) stop("readLength must be specified")
  if (class(genomeDB)!='annotatedGenome') stop("genomeDB must be of class 'annotatedGenome'")
  if (!genomeDB@denovo) stop("genomeDB must be a de novo annotated genome. Use createDenovoGenome")
  if (class(pc)!="pathCounts") stop("pc must be of class 'pathCounts'")
  if (!pc@denovo) stop("pc must be computed using a genome annotated de novo")
  if (!missing(mprior)) {
    if (!all(c('nvarPrior','nexonPrior') %in% slotNames(mprior))) stop("Incorrect mprior. Please use modelPrior to generate it.")
    modelUnifPrior <- as.integer(0)
    nvarPrior <- as.list(data.frame(t(mprior@nvarPrior$nbpar)))
    nexonPrior <- as.list(data.frame(t(mprior@nexonPrior$bbpar)))
  } else {
    nvarPrior <- list(nbpar=matrix(c(0,0),nrow=1),obs=NA,pred=NA)
    nexonPrior <- list(bbpar=matrix(c(0,0),nrow=1),obs=NA,pred=NA)
    modelUnifPrior <- as.integer(1)
  }
  if (!(method %in% c('auto','rwmcmc','priormcmc','exact'))) stop("method must be auto, rwmcmc, priormcmc or exact")
  
  #Format input
  startcdf <- as.double(ecdf(distrs$stDis)(seq(0,1,.001)))
  lendis <- as.double(distrs$lenDis/sum(distrs$lenDis))
  lenvals <- as.integer(names(distrs$lenDis))
  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  nvarPrior <- lapply(nvarPrior,as.double)
  nexonPrior <- lapply(nexonPrior,as.double)
  minpp <- as.double(minpp)
  selectBest <- as.integer(selectBest)
  if (method=='auto') method <- 0 else if (method=='exact') method <- 1 else if (method=='rwmcmc') method <- 2 else method <- 3
  method <- as.integer(method)
  verbose <- as.integer(verbose)
  niter <- as.integer(niter)
  exactMarginal <- as.integer(exactMarginal)
  if (missing(geneid)) geneid <- names(genomeDB@islands)[sapply(genomeDB@islands,length)>1]

  exons <- lapply(genomeDB@islands,function(z) as.integer(names(z)))
  exonwidth <- lapply(genomeDB@islands,width)

  if (!all(geneid %in% names(exons))) stop('geneid not found in genomeDB@islands')
  if (!all(geneid %in% names(pc@counts))) stop('geneid not found in pc')
  if (!all(geneid %in% names(genomeDB@transcripts))) stop('geneid not found in genomeDB@transcripts')

  #Define basic function
  f <- function(z) {
    geneid <- as.integer(z)
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- genomeDB@transcripts[z]
    pc <- pc@counts[z]
    ans <- calcDenovoMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,geneid=as.list(geneid),pc=pc,startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,modelUnifPrior=modelUnifPrior,nvarPrior=nvarPrior,nexonPrior=nexonPrior,priorq=priorq,minpp=minpp,selectBest=selectBest,method=method,niter=niter,exactMarginal=exactMarginal,verbose=verbose)
    mapply(function(z1,z2) formatDenovoOut(z1,z2), ans, genomeDB@islands[z], SIMPLIFY=FALSE)
  }

  runCalc <- function(geneid) {
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
    ans
  }


  #Initialize transcripts for new islands with known orientation
  sel <- names(genomeDB@transcripts)[sapply(genomeDB@transcripts,is.null) & !is.na(genomeDB@txStrand)]
  if (length(sel)>0) genomeDB@transcripts[sel] <- lapply(genomeDB@exons[sel],function(z) list(as.numeric(names(z))))

  geneidUnknown <- geneid[geneid %in% names(genomeDB@transcripts)[sapply(genomeDB@transcripts,is.null)]]
  if (length(geneidUnknown)>0) { geneidini <- geneid; geneid <- geneid[!(geneid %in% geneidUnknown)] }

  #Run
  if (length(geneidUnknown)==0) {
    ans <- runCalc(geneid)
  } else {
    #Islands with known strand
    ans <- vector("list",length(geneidini)); names(ans) <- geneidini
    if (length(geneid)>0) ans[geneid] <- runCalc(geneid)
    
    #Islands with unknown strand. Run 2 strands and select the one with largest post prob
    genomeDB@transcripts[geneidUnknown] <- lapply(genomeDB@islands[geneidUnknown],function(z) list(var1=as.integer(names(z))))
    ansforw <- runCalc(geneidUnknown)
    genomeDB@transcripts[geneidUnknown] <- lapply(genomeDB@islands[geneidUnknown],function(z) list(var1=rev(as.integer(names(z)))))
    ansrev <- runCalc(geneidUnknown)
    difndel <- sapply(ansforw,function(z) z@npathDeleted) - sapply(ansrev,function(z) z@npathDeleted)
    difmax <- sapply(ansforw,function(z) z@integralSum['logmax']) - sapply(ansrev,function(z) z@integralSum['logmax'])
    difsum <- log(sapply(ansforw,function(z) z@integralSum['sum'])) - log(sapply(ansrev,function(z) z@integralSum['sum']))
    sel <- ifelse(difndel<0 || (difndel==0 && (difsum+difmax)>=0), TRUE, FALSE)
    ans[geneidUnknown[sel]] <- ansforw[sel]; ans[geneidUnknown[!sel]] <- ansrev[!sel]
  }
  
  new("denovoGenomeExpr", islands=ans)
}



formatDenovoOut <- function(ans, genesel) {
  ans[[1]] <- data.frame(ans[[1]],ans[[6]])
  ans[[6]] <- NULL
  colnames(ans[[1]]) <- c('model','posprob','modelid')
  ans[[2]] <- data.frame(ans[[2]],ans[[3]])
  ans[[3]] <- NULL
  colnames(ans[[2]]) <- c('model','expr','varName')
  variants <- genesel[as.character(ans[[3]])]
  names(variants) <- NULL
  ans[[3]] <- RangedData(ranges=variants,space=ans[[4]])
  ans[[4]] <- NULL
  names(ans[[4]]) <- c('sum','logmax')
  names(ans) <- c('posprob','expression','variants','integralSum','npathDeleted')
  new("denovoGeneExpr",posprob=ans$posprob,expression=ans$expression,variants=ans$variants,integralSum=ans$integralSum,npathDeleted=ans$npathDeleted)
}


calcDenovoMultiple <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, modelUnifPrior, nvarPrior, nexonPrior, priorq, minpp, selectBest, method, niter, exactMarginal, verbose) {
  ans <- .Call("calcDenovoMultiple",exons,exonwidth,transcripts,geneid,pc,startcdf,lendis,lenvals,readLength,modelUnifPrior,nvarPrior,nexonPrior,priorq,minpp,selectBest,method,niter,exactMarginal,verbose)
  return(ans)
}


