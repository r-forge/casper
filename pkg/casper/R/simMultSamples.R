simMultSamples <- function(B, nsamples, nreads, readLength, x, groups='group', distrs, genomeDB, verbose=TRUE, mc.cores=1) {
# Posterior predictive simulation for multiple samples
# nsamples: vector w/ n. samples per group
# nreads: nreads per sample
# readLength: read length
# x: ExpressionSet pilot data. x[[group]] indicates groups to be compared
# groups: name of column in pData(x) indicating the groups
# distrs: fragment start and length distributions. It can be either an object of class readDistrs, or a list where each element is of class readDistrs. In the latter case, an element is chosen at random for each individual sample (so that uncertainty in these distributions is taken into account).
# genomeDB: annotatedGenome object 
  require(gaga)
  require(plyr)
  if (verbose) cat("Fitting NNGV model...\n")
  seed <- sample(1:10000, 1)
  l <- genomeDB@txLength
  groupsnew <- rep(unique(pData(x)[,groups]), nsamples)
  nnfit <- fitNNSingleHyp(x, groups=groups, B=5, trace=FALSE)
  ans <- vector("list",B)
  if (verbose) cat(paste("Obtaining ",B," simulations (",sum(nsamples)," samples with ",nreads," reads each)\n",sep=''))
  for (k in 1:B) {
    xnew <- simnewsamples(nnfit, groupsnew=groupsnew, x=x, groups=groups)
    sampleNames(xnew) <- paste("Sample",1:ncol(xnew))
    rownames(fData(xnew)) <- featureNames(xnew)
    # fData(xnew), simulated (phi, mu1, mu2)
    # xnew: simulated (unobservable) expressions for new individuals 
    rownames(exprs(xnew)) <- rownames(exprs(x)) 
    isof.xnew <- rownames(exprs(xnew))                
    included <- isof.xnew[isof.xnew %in% names(l)] # this step might be innecessary  
    thpi <- exp(exprs(xnew)[included,])*l[included]
    thpi <- t(t(thpi)/colSums(thpi))
    g <- as.numeric(getIsland(txid=rownames(thpi), genomeDB=genomeDB))
    th <- aggregate(thpi, by=list(g), FUN=sum)[,-1]
    th2 <- as.data.frame(th)
    th2$islandid <- as.numeric(rownames(th))
    th.ext <- data.frame(islandid=g)
    th2 <- join(th.ext, th2, by='islandid')
    th2 <- as.matrix(th2)
    rownames(th2) <- th2[,1]
    th2 <-  th2[,-1]
    pis <- thpi/th2
    N <- apply(th, 2, function(x) rmultinom(size=nreads, n=1, prob=x))
    rownames(N) <- rownames(th)
    if (mc.cores>1) {
      if ('parallel' %in% loadedNamespaces()) {
        sim.exp <- parallel::mclapply(1:ncol(xnew), simOneExp, mc.cores=mc.cores, mc.preschedule=FALSE)
      } else {
        stop("parallel package has not been loaded!")
      }
    } else {
      sim.exp <- lapply(1:ncol(xnew), simOneExp)
    }
    explCnts <- do.call(cbind,lapply(sim.exp,'[[','explCnts'))
    colnames(explCnts) <- paste('explCnts',1:ncol(explCnts))
    sim.exp <- do.call(cbind,lapply(sim.exp,'[[','exp'))
    colnames(sim.exp) <- sampleNames(xnew)
    b <- new("AnnotatedDataFrame", data=pData(xnew))
    e <- new("ExpressionSet", exprs=sim.exp, phenoData=b, featureData=new("AnnotatedDataFrame", explCnts)) 
    ans[[i]] <- list(simTruth=fData(xnew), simExpr=e)
    if (verbose) cat("\n")
  }
  return(ans)
} 


simOneExp <- function(i) {
  #Simulate a single sample
  if (is.list(distrs)) distrsCur <- sample(distrs, 1)[[1]] else distrsCur <- distrs
  nSimReads <- N[which(N[,i] != 0),i]
  islandid <- names(nSimReads)
  zero <- which(N[,i]==0)
  sim.r <- simReads(writeBam=0, seed=seed, islandid=islandid, nSimReads=nSimReads, pis=pis[,i], rl=as.integer(readLength), distrs=distrsCur, genomeDB=genomeDB, verbose=FALSE)
  v.list <- list(NULL)
  length(v.list) <- length(zero)
  names(v.list) <-  as.character(zero)
  sim.r@counts[[1]] <- c(v.list, sim.r@counts[[1]])
  exp.sim <- calcExp(distrs=distrsCur, genomeDB=genomeDB, pc=sim.r, readLength=readLength)
  ans <- list(explCnts=fData(exp.sim)[featureNames(xnew),'explCnts',drop=FALSE], exp=exprs(exp.sim)[featureNames(xnew),,drop=FALSE])
  if (verbose) cat('.')
  return(ans)
}
