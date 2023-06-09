#########################################################################
## DEFINITION AND METHODS FOR CLASS denovoGeneExpr AND denovoGenomeExpr
#########################################################################

require(methods)

setClass("denovoGeneExpr", representation(posprob = "data.frame", expression = "data.frame", variants = "GRanges", integralSum= "numeric", npathDeleted= "numeric"))

valid_denovoGeneExpr <- function(object) {
  msg <- NULL
  if (any(!(c('model','posprob') %in% colnames(object@posprob)))) msg <- "Incorrect colnames in 'posprob'"
  if (any(!(c('model','expr','varName') %in% colnames(object@expression)))) msg <- "Incorrect colnames in 'expression'"
  if (class(object@variants)!='GRanges') msg <- "Element variants must be of class 'Granges'"
  if (is.null(msg)) { TRUE } else { msg }
}

setValidity("denovoGeneExpr", valid_denovoGeneExpr)

setMethod("show", signature(object="denovoGeneExpr"), function(object) {
  cat("denovoGeneExpr object\n")
  cat("\nPosterior model probabilities\n")
  show(head(object@posprob))
  cat("...\nEstimated expression (conditional on each model)\n")
  show(head(object@expression))
  cat("...\nUse posprob() to access posterior probabilities; variants() to get exons in each variant\n")
}
)

setGeneric("posprob", function(object) standardGeneric("posprob"))
setMethod("posprob", signature(object="denovoGeneExpr"), function(object) {
  object@posprob
}
)

setGeneric("variants", function(object) standardGeneric("variants"))
setMethod("variants", signature(object="denovoGeneExpr"), function(object) {
  object@variants
}
)

setGeneric("variants<-", function(object,value) standardGeneric("variants<-"))
setReplaceMethod("variants", "denovoGeneExpr", function(object, value) { object@variants <- value; object })



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
setMethod("as.list", signature(x="denovoGenomeExpr"), function(x) {x@islands})
                            

#########################################################################
## Function calcDenovo
#########################################################################

calcDenovo <- function(distrs, genomeDB, pc, readLength, islandid, priorq=3, mprior, minpp=0.001, selectBest=FALSE, method='auto', niter, exactMarginal=TRUE, verbose=TRUE, mc.cores=1) {
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
  if (verbose) cat("Formatting input...\n")
  sseq <- seq(0,1,.001)
  startcdf <- as.double(distrs@stDis(sseq))

  lenvals <- as.integer(names(distrs@lenDis))
  lenvals <- as.integer(seq(min(lenvals),max(lenvals),1))
  lendis <- rep(0,length(lenvals)); names(lendis) <- as.character(lenvals)
  lendis[names(distrs@lenDis)] <- as.double(distrs@lenDis/sum(distrs@lenDis))

  readLength <- as.integer(readLength)
  priorq <- as.double(priorq)
  nvarPrior <- lapply(nvarPrior,as.double)
  nexonPrior <- lapply(nexonPrior,as.double)
  minpp <- as.double(minpp)
  selectBest <- as.integer(selectBest)
  if (method=='auto') method <- 0 else if (method=='exact') method <- 1 else if (method=='rwmcmc') method <- 2 else method <- 3
  method <- as.integer(method)
  verbose <- as.integer(verbose)
  exactMarginal <- as.integer(exactMarginal)
  if (missing(islandid)) {
    islandid <- names(genomeDB@islands)[elementLengths(genomeDB@islands)>1]
    islandid <- islandid[islandid %in% names(pc@counts[[1]])]
  }
  if (!all(islandid %in% names(genomeDB@islands))) stop('islandid not found in genomeDB@islands')
  if (!all(islandid %in% unlist(lapply(pc@counts, names)))) stop('islandid not found in pc')
  if (!all(islandid %in% names(genomeDB@transcripts))) stop('islandid not found in genomeDB@transcripts')
  exons <- as.integer(names(genomeDB@islands@unlistData))
  names(exons) <- rep(names(genomeDB@islands), elementLengths(genomeDB@islands))
  exons <- split(unname(exons), names(exons))
  exonwidth <- width(genomeDB@islands@unlistData)
  names(exonwidth) <- rep(names(genomeDB@islands), elementLengths(genomeDB@islands))
  exonwidth <- split(unname(exonwidth), names(exonwidth))
  strand <- as.character(strand(genomeDB@islands@unlistData))[cumsum(c(1, elementLengths(genomeDB@islands)[-length(genomeDB@islands)]))]
  names(strand) <- names(genomeDB@islands)
  
  if (missing(niter)) {
     niter <- as.list(as.integer(ifelse(sapply(exons,length)>20,10^3,10^4)))
  } else {
     niter <- as.list(as.integer(rep(niter[1],length(islandid))))
  }
  names(niter) <- islandid
  
  #Define basic function
  f <- function(z) {
    islandid <- as.integer(as.factor(z))
    exons <- exons[z]
    exonwidth <- exonwidth[z]
    transcripts <- genomeDB@transcripts[z]
    tmp <- strand[z]
    strand <- vector(mode='integer', length=length(tmp))
    sel <- tmp=='+'
    strand[sel] <- 1
    sel <- tmp=='-'
    strand[sel] <- -1
    sel <- tmp=='*'
    strand[sel] <- 0
    strand <- as.list(as.integer(strand))
    #pc <- pc[z] pc's have to be subset from previous step to deal with strandedness !!!!!!
    ans <- calcDenovoMultiple(exons=exons,exonwidth=exonwidth,transcripts=transcripts,islandid=as.list(islandid),pc=pc@counts[[1]][z],startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,modelUnifPrior=modelUnifPrior,nvarPrior=nvarPrior,nexonPrior=nexonPrior,priorq=priorq,minpp=minpp,selectBest=selectBest,method=method,niter=niter[z],exactMarginal=exactMarginal,verbose=verbose, strand=strand)
    lapply(1:length(ans), function(y) formatDenovoOut(ans[[y]], genomeDB@islands[z][[y]]))
  }

  formatZeroExpr <- function(ids){
    isl <- genomeDB@islands[ids]
    txs <- genomeDB@transcripts[ids]
    spa <- lapply(txs, function(x) rep(names(x), unlist(lapply(x, length))))
      #exo <- lapply(names(txs), function(x) {  tmp <- genomeDB@islands[[x]][as.character(unlist(txs[[x]]))]; names(tmp) <- sub("\\..*", "", names(tmp)); tmp})
    #names(exo) <- names(txs)
    exo <- isl
    expr <- lapply(ids, function(x) data.frame(model=rep(0, length(exo[[x]])), expr=rep(0, length(exo[[x]])), varName=names(exo[[x]])))
    names(expr) <- ids
    posprob <- lapply(ids, function(x) data.frame(model=0, posprob=NA, modelid=0))
    names(posprob) <- ids
    res <- lapply(ids, function(x) new("denovoGeneExpr", variants=exo[[x]], expression=expr[[x]], posprob=posprob[[x]]))
    names(res) <- ids
    res
  }
    
  runCalc <- function(islandid) {
    sel <- !(sapply(pc@counts[[1]][islandid], is.null) | strand[islandid]=='*')
    if (verbose==1) cat("Note: Islands with transcripts from both strands will not be processed at the moment\n")
    all <- islandid
    islandid <- islandid[sel]
    if(sum(sel)==0) stop("Not islands left to process due to strand or lack of path counts\n")
    if (mc.cores>1 && length(islandid)>mc.cores) {
      require(parallel)
      if ('parallel' %in% loadedNamespaces()) {
        #ans <- mclapply(islandid, f, mc.cores=mc.cores)
        #split into smaller jobs
        nsplit <- ceiling(max(length(islandid), mc.cores)/mc.cores)
        islandidList <- lapply(1:min(length(islandid), mc.cores), function(z) islandid[seq(z,length(islandid),by=mc.cores)])
        ans <- parallel::mclapply(islandidList,f,mc.cores=min(length(islandidList), mc.cores))
        ans <- do.call(c,ans); names(ans) <- unlist(islandidList); ans <- ans[islandid]
      } else stop('parallel library has not been loaded!')
    } else {
      ans <- f(islandid)
      names(ans) <- islandid
    }
    res <- vector(mode="list", length=length(all))
    names(res) <- all
    z <- formatZeroExpr(all[!sel])
    res[names(z)] <- z
    res[names(ans)] <- ans
    res
  }

  #Initialize transcripts for new islands with known orientation
  sel <- names(genomeDB@transcripts)[sapply(genomeDB@transcripts,is.null) & strand=='*']
  if (length(sel)>0) genomeDB@transcripts[sel] <- tapply(as.integer(names(genomeDB@islands[sel]@unlistData)), rep(names(genomeDB@islands[sel]), elementLengths(genomeDB@islands[sel])), function(x) list(as.numeric(x))) 
  islandidUnknown <- islandid[islandid %in% names(genomeDB@transcripts)[sapply(genomeDB@transcripts,is.null)]]
  if (length(islandidUnknown)>0) { islandidini <- islandid; islandid <- islandid[!(islandid %in% islandidUnknown)] }


  #Run
  if (verbose==1) cat("Performing model search (this may take a while)\n")
  if (length(islandidUnknown)==0) {
    ans <- runCalc(islandid)
  } else {
    #Islands with known strand
    ans <- vector("list",length(islandidini)); names(ans) <- islandidini
    if (length(islandid)>0) ans[islandid] <- runCalc(islandid)
    
    #Islands with unknown strand. Run 2 strands and select the one with largest post prob
    genomeDB@transcripts[islandidUnknown] <- lapply(genomeDB@islands[islandidUnknown],function(z) list(var1=as.integer(names(z))))
    strand[islandidUnknown] <- '+'
    ansforw <- runCalc(islandidUnknown)
    strand[islandidUnknown] <- '-'
    genomeDB@transcripts[islandidUnknown] <- lapply(genomeDB@transcripts[islandidUnknown],rev)
    genomeDB@islands[islandidUnknown] <- lapply(genomeDB@islands[islandidUnknown], rev)
    exons[islandidUnknown] <- lapply(exons[islandidUnknown], rev)
    exonwidth[islandidUnknown] <- lapply(exonwidth[islandidUnknown], rev)
    ansrev <- runCalc(islandidUnknown)
    difndel <- sapply(ansforw,function(z) z@npathDeleted) - sapply(ansrev,function(z) z@npathDeleted)
    difmax <- sapply(ansforw,function(z) z@integralSum['logmax']) - sapply(ansrev,function(z) z@integralSum['logmax'])
    difsum <- log(sapply(ansforw,function(z) z@integralSum['sum'])) - log(sapply(ansrev,function(z) z@integralSum['sum']))
    sel <- ifelse(difndel<0 || (difndel==0 && (difsum+difmax)>=0), TRUE, FALSE)
    ans[islandidUnknown[sel]] <- ansforw[sel]; ans[islandidUnknown[!sel]] <- ansrev[!sel]
  }
  if (verbose==1) cat("\n")
  new("denovoGenomeExpr", islands=ans)
}


formatDenovoOut <- function(ans, genesel) {
  ans[[1]] <- data.frame(ans[[1]])
  colnames(ans[[1]]) <- c('model','posprob','priorprob')
  ans[[1]] <- ans[[1]][order(ans[[1]][,'posprob'],decreasing=TRUE),]
  ans[[6]] <- NULL
  ans[[2]] <- data.frame(ans[[2]],ans[[3]])
  colnames(ans[[2]]) <- c('model','expr','varName')
  ans[[3]] <- NULL
  variants <- genesel[as.character(ans[[3]])]
  names(variants) <- sub("\\..*", "", names(variants))
  ans[[3]] <- variants
  ans[[4]] <- NULL
  names(ans[[4]]) <- c('sum','logmax')
  names(ans) <- c('posprob','expression','variants','integralSum','npathDeleted')
  new("denovoGeneExpr",posprob=ans$posprob,expression=ans$expression,variants=ans$variants,integralSum=ans$integralSum,npathDeleted=ans$npathDeleted)
}

calcDenovoMultiple <- function(exons, exonwidth, transcripts, islandid, pc, startcdf, lendis, lenvals, readLength, modelUnifPrior, nvarPrior, nexonPrior, priorq, minpp, selectBest, method, niter, exactMarginal, verbose, strand) {
  ans <- .Call("calcDenovoMultiple",exons,exonwidth,transcripts,islandid,pc,startcdf,lendis,lenvals,readLength,modelUnifPrior,nvarPrior,nexonPrior,priorq,minpp,selectBest,method,niter,exactMarginal,verbose, strand)
  return(ans)
}



variantMargExpr <- function(x,minProbExpr=0.5, minExpr=0.05) {
  #Marginal expression for each variant (obtained via model averaging) and marginal post prob of being expressed
  # - minProbExpr: variants with marginal post prob < minProbExpr are not reported
  # - minExpr: variants with expression < minExpr are not reported
  # Note: at least one variant is always reported. If no variants satisfy minProbExpr and minExpr, the variant with largest expression is reported
  pospr <- x@posprob$posprob/sum(x@posprob$posprob)
  names(pospr) <- x@posprob$model
  pospr <- pospr[as.character(x@expression$model)]
  ans <- by(data.frame(pospr*x@expression$expr,pospr),INDICES=list(var=x@expression$varName),FUN=colSums,simplify=FALSE)
  n <- names(ans)
  ans <- matrix(unlist(ans),ncol=2,byrow=TRUE)
  colnames(ans) <- c('expr','probExpressed')
  rownames(ans) <- n
  sel <- ans[,'probExpressed']>minProbExpr & ans[,'expr']>minExpr
  if (any(sel)) ans <- ans[sel,,drop=FALSE] else ans <- ans[which.max(ans[,'expr']),,drop=FALSE]
  ans[,'expr'] <- ans[,'expr']/sum(ans[,'expr'])
  return(ans)
}

relativeExpr <- function(expr, summarize='modelAvg', minProbExpr=0.5, minExpr=0.05){
  if (class(expr)!='denovoGenomeExpr') stop("expr must be of class 'denovoGenomeExpr'")
  if (!(summarize %in% c("bestModel", "modelAvg"))) stop("summarize must be one of 'bestModel' or 'modelAvg'")
  if (summarize=='bestModel'){
    ans <- lapply(as.list(expr), function(x){
      best <- x@posprob$model[which.max(x@posprob$posprob)]
      exp <- x@expression[x@expression$model==best,]
      res <- exp$expr
      names(res) <- exp$varName
      res
    })
    ans <- do.call(c, unname(ans))
  } else {
    ans <- lapply(as.list(expr), variantMargExpr, minProbExpr=minProbExpr, minExpr=minExpr)
    ans <- do.call("rbind", unname(ans))
  }
  ans
}


