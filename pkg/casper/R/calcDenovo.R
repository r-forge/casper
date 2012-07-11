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

calcDenovo <- function(distrs, genomeDB, pc, readLength, geneid, priorq=3, mprior, minpp=0.001, selectBest=FALSE, method='auto', niter, exactMarginal=TRUE, verbose=TRUE, mc.cores=1) {
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
  startcdf <- as.double(distrs$stDis(sseq))

  lenvals <- as.integer(names(distrs$lenDis))
  lenvals <- as.integer(seq(min(lenvals),max(lenvals),1))
  lendis <- rep(0,length(lenvals)); names(lendis) <- as.character(lenvals)
  lendis[names(distrs$lenDis)] <- as.double(distrs$lenDis/sum(distrs$lenDis))

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
  if (missing(geneid)) geneid <- names(genomeDB@islands)[sapply(genomeDB@islands,length)>1]

  if (!all(geneid %in% names(genomeDB@islands))) stop('geneid not found in genomeDB@islands')
  if (!all(geneid %in% names(pc@counts))) stop('geneid not found in pc')
  if (!all(geneid %in% names(genomeDB@transcripts))) stop('geneid not found in genomeDB@transcripts')
  
  exons <- lapply(genomeDB@islands[geneid],function(z) as.integer(names(z)))
  exonwidth <- lapply(genomeDB@islands[geneid],width)
  strand <- genomeDB@islandStrand
  if (missing(niter)) {
     niter <- as.list(as.integer(ifelse(sapply(exons,length)>20,10^3,10^4)))
  } else {
     niter <- as.list(as.integer(rep(niter[1],length(geneid))))
  }
  names(niter) <- geneid
  
  #Define basic function
  f <- function(z) {
    geneid <- as.integer(z)
    transcripts <- genomeDB@transcripts[z]
    strand <- as.list(as.integer(ifelse(strand[z]=='+', 1,-1)))
    ans <- calcDenovoMultiple(exons=exons[z],exonwidth=exonwidth[z],transcripts=transcripts,geneid=as.list(geneid),pc=pc@counts[z],startcdf=startcdf,lendis=lendis,lenvals=lenvals,readLength=readLength,modelUnifPrior=modelUnifPrior,nvarPrior=nvarPrior,nexonPrior=nexonPrior,priorq=priorq,minpp=minpp,selectBest=selectBest,method=method,niter=niter[z],exactMarginal=exactMarginal,verbose=verbose, strand=strand)
    mapply(function(z1,z2) formatDenovoOut(z1,z2), ans, genomeDB@islands[z], SIMPLIFY=FALSE)
  }

  formatZeroExpr <- function(ids){
    isl <- genomeDB@islands[ids]
    txs <- genomeDB@transcripts[ids]
    spa <- lapply(txs, function(x) rep(names(x), unlist(lapply(x, length))))
    exo <- lapply(names(txs), function(x) genomeDB@islands[[x]][as.character(unlist(txs[[x]])),])
    names(exo) <- names(txs)
    ran <- lapply(names(txs), function(x) RangedData(unname(exo[[x]]), space=unlist(spa[[x]])))
    names(ran) <- ids
    expr <- lapply(ids, function(x) data.frame(model=rep(0, length(ran[[x]])), expr=rep(0, length(ran[[x]])), varName=names(ran[[x]])))
    names(expr) <- ids
    posprob <- lapply(ids, function(x) data.frame(model=0, posprob=NA, modelid=0))
    names(posprob) <- ids
    res <- lapply(ids, function(x) new("denovoGeneExpr", variants=ran[[x]], expression=expr[[x]], posprob=posprob[[x]]))
    names(res) <- ids
    res
  }
    
  runCalc <- function(geneid) {
    sel <- !sapply(pc@counts[geneid], is.null)
    all <- geneid
    geneid <- geneid[sel]
    if (mc.cores>1 && length(geneid)>mc.cores) {
      if ('multicore' %in% loadedNamespaces()) {
        #ans <- mclapply(geneid, f, mc.cores=mc.cores)
        #split into smaller jobs
        nsplit <- ceiling(max(length(geneid), mc.cores)/mc.cores)
        geneidList <- lapply(1:min(length(geneid), mc.cores), function(z) geneid[seq(z,length(geneid),by=mc.cores)])
        ans <- mclapply(geneidList,f,mc.cores=min(length(geneidList), mc.cores))
        ans <- do.call(c,ans); names(ans) <- unlist(geneidList); ans <- ans[geneid]
      } else stop('multicore library has not been loaded!')
    } else {

      ans <- f(geneid)
      names(ans) <- geneid
    }
    res <- vector(mode="list", length=length(all))
    names(res) <- all
    z <- formatZeroExpr(all[!sel])
    res[names(z)] <- z
    res[names(ans)] <- ans
    res
  }


  #Initialize transcripts for new islands with known orientation
  sel <- names(genomeDB@transcripts)[sapply(genomeDB@transcripts,is.null) & !is.na(genomeDB@islandStrand)]
  if (length(sel)>0) genomeDB@transcripts[sel] <- lapply(genomeDB@exons[sel],function(z) list(as.numeric(names(z))))

  geneidUnknown <- geneid[geneid %in% names(genomeDB@transcripts)[sapply(genomeDB@transcripts,is.null)]]
  if (length(geneidUnknown)>0) { geneidini <- geneid; geneid <- geneid[!(geneid %in% geneidUnknown)] }

  #Run
  if (verbose==1) cat("Performing model search (this may take a while)")
  if (length(geneidUnknown)==0) {
    ans <- runCalc(geneid)
  } else {
    #Islands with known strand
    ans <- vector("list",length(geneidini)); names(ans) <- geneidini
    if (length(geneid)>0) ans[geneid] <- runCalc(geneid)
    
    #Islands with unknown strand. Run 2 strands and select the one with largest post prob
    genomeDB@transcripts[geneidUnknown] <- lapply(genomeDB@islands[geneidUnknown],function(z) list(var1=as.integer(names(z))))
    strand[geneidUnknown] <- '+'
    ansforw <- runCalc(geneidUnknown)
    strand[geneidUnknown] <- '-'
    genomeDB@transcripts[geneidUnknown] <- lapply(genomeDB@transcripts[geneidUnknown],rev)
    genomeDB@islands[geneidUnknown] <- lapply(genomeDB@islands[geneidUnknown], rev)
    exons[geneidUnknown] <- lapply(exons[geneidUnknown], rev)
    exonwidth[geneidUnknown] <- lapply(exonwidth[geneidUnknown], rev)
    ansrev <- runCalc(geneidUnknown)
    difndel <- sapply(ansforw,function(z) z@npathDeleted) - sapply(ansrev,function(z) z@npathDeleted)
    difmax <- sapply(ansforw,function(z) z@integralSum['logmax']) - sapply(ansrev,function(z) z@integralSum['logmax'])
    difsum <- log(sapply(ansforw,function(z) z@integralSum['sum'])) - log(sapply(ansrev,function(z) z@integralSum['sum']))
    sel <- ifelse(difndel<0 || (difndel==0 && (difsum+difmax)>=0), TRUE, FALSE)
    ans[geneidUnknown[sel]] <- ansforw[sel]; ans[geneidUnknown[!sel]] <- ansrev[!sel]
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

calcDenovoMultiple <- function(exons, exonwidth, transcripts, geneid, pc, startcdf, lendis, lenvals, readLength, modelUnifPrior, nvarPrior, nexonPrior, priorq, minpp, selectBest, method, niter, exactMarginal, verbose, strand) {
  ans <- .Call("calcDenovoMultiple",exons,exonwidth,transcripts,geneid,pc,startcdf,lendis,lenvals,readLength,modelUnifPrior,nvarPrior,nexonPrior,priorq,minpp,selectBest,method,niter,exactMarginal,verbose, strand)
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

relativeExpr <- function(expr, method='modelAvg', minProbExpr=0.5, minExpr=0.05){
  if (!(method %in% c("bestModel", "modelAvg"))) stop("method must be one of 'bestModel' or 'modelAvg'")
  if (class(expr)!='denovoGenomeExpr') stop("expr must be of class 'denovoGenomeExpr'")
  if (method=='bestModel'){
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
    ans <- do.call("c", unname(ans))
  }
  ans
}


