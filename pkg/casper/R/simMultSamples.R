logrpkm2thpi <- function(txids, lrpkm, genomeDB) {
  #txids: data.frame with transcript, gene
  #lrpkm: matrix with expression levels measured in logrpkm. rownames must indicate tx name
  #genomeDB: annotatedGenome
  rownames(txids) <- as.character(txids$transcript)
  txids <- txids[rownames(lrpkm),]
  len <- genomeDB@txLength[rownames(lrpkm)]
  thpi <- data.frame(txids[,c('transcript','gene')], exp(lrpkm) * len / 10^9)
  thpi[,-1:-2] <- t(t(thpi[,-1:-2,drop=FALSE])/colSums(thpi[,-1:-2,drop=FALSE]))
  th <- aggregate(thpi[,-1:-ncol(txids),drop=FALSE], by=list(thpi$gene), FUN=sum)
  names(th)[1] <- 'gene'
  thpi <- merge(thpi, th, by='gene')
  pi <- thpi[,3:(2+ncol(lrpkm)),drop=FALSE] / thpi[,(3+ncol(lrpkm)):ncol(thpi),drop=FALSE]
  rownames(pi) <- thpi$transcript; colnames(pi) <- colnames(lrpkm)
  pi <- pi[rownames(lrpkm),]
  return(list(th=th,pi=pi))
}


logrpkm <- function(txids, th, pi, genomeDB) {
  #txids: data.frame with transcript, gene
  #th: data.frame with gene and proportion/probability of reads from each gene
  #pi: matrix with relative expression of each isoform, rownames must indicate tx identifier
  #genomeDB: annotatedGenome
  txids$gene <- as.character(txids$gene)
  th$gene <- as.character(th$gene)
  th <- merge(txids[,c('transcript','gene')], th, by='gene')
  rownames(th) <- as.character(th$transcript)
  th <- th[rownames(pi),,drop=FALSE]
  len <- genomeDB@txLength[rownames(pi)]
  lrpkm <- log(10^9) + log(th[,-1:-2,drop=FALSE]) + log(pi) - log(len)
  return(lrpkm)
}


simMultSamples <- function(B, nsamples, nreads, readLength, x, groups='group', distrs, genomeDB, verbose=TRUE, mc.cores=1) {
# Posterior predictive simulation for multiple samples
# nsamples: vector w/ n. samples per group
# nreads: nreads per sample
# readLength: read length
# x: ExpressionSet pilot data. x[[group]] indicates groups to be compared
# groups: name of column in pData(x) indicating the groups
# distrs: fragment start and length distributions. It can be either an object of class readDistrs, or a list where each element is of class readDistrs. In the latter case, an element is chosen at random for each individual sample (so that uncertainty in these distributions is taken into account).
# genomeDB: annotatedGenome object
  if (!all(paste(sampleNames(x),'se',sep='.') %in% names(fData(x)))) {
    stop("Estimate SE cannot be found in fData(x). Please run calcExp with citype='asymp'")
  } else {
    sigma2ErrorObs <- fData(x)[,paste(sampleNames(x),'se',sep='.'),drop=FALSE]^2
  }
  if (verbose) cat("Fitting NNGV model...\n")
  seed <- sample(1:10000, 1)
  l <- genomeDB@txLength[featureNames(x)]
  groupsnew <- rep(unique(pData(x)[,groups]), nsamples)
  nnfit <- fitNNSingleHyp(x, groups=groups, B=5, trace=FALSE)
  ans <- vector("list",B)
  if (verbose) cat(paste("Obtaining ",B," simulations (",sum(nsamples)," samples with ",nreads," reads each)\n",sep=''))
  for (k in 1:B) {
    xnew <- simnewsamplesNoisyObs(nnfit, groupsnew=groupsnew, x=x, groups=groups, sigma2ErrorObs=sigma2ErrorObs)
    #xnew <- t(t(exprs(xnew)) - colMeans(exprs(xnew)) + colMeans(exprs(x)))
    sampleNames(xnew) <- paste("Sample",1:ncol(xnew))
    featureNames(xnew) <- featureNames(x)
    # fData(xnew), simulated (phi, mu1, mu2)
    # xnew: simulated (unobservable) expressions for new individuals
    thpi <- logrpkm2thpi(txids=fData(x)[,c('transcript','gene')], lrpkm=exprs(xnew), genomeDB=genomeDB)
    N <- apply(thpi$th[,-1], 2, function(x) rmultinom(size=nreads, n=1, prob=x))
    rownames(N) <- as.character(thpi$th$gene)
    if (mc.cores>1) {
      if ('parallel' %in% loadedNamespaces()) {
        sim.exp <- parallel::mclapply(1:ncol(xnew), simOneExp, distrs=distrs, N=N, pis=thpi$pi, readLength=readLength, genomeDB=genomeDB, featnames=featureNames(xnew), seed=seed, verbose=verbose, mc.cores=mc.cores, mc.preschedule=FALSE)
      } else {
        stop("parallel package has not been loaded!")
      }
    } else {
      sim.exp <- lapply(1:ncol(xnew), simOneExp, distrs=distrs, N=N, pis=pis, readLength=readLength, genomeDB=genomeDB, featnames=featureNames(xnew), seed=seed)
    }
    explCnts <- do.call(cbind,lapply(sim.exp,'[[','explCnts'))
    colnames(explCnts) <- paste('explCnts',1:ncol(explCnts))
    sim.exp <- do.call(cbind,lapply(sim.exp,'[[','exp'))
    colnames(sim.exp) <- sampleNames(xnew)
    b <- new("AnnotatedDataFrame", data=pData(xnew))
    e <- new("ExpressionSet", exprs=sim.exp, phenoData=b, featureData=new("AnnotatedDataFrame", explCnts)) 
    ans[[k]] <- list(simTruth=fData(xnew), simExprTrue=exprs(xnew), simExpr=e)
    if (verbose) cat("\n")
  }
  return(ans)
} 


simOneExp <- function(i, distrs, N, pis, readLength, genomeDB, featnames, seed, verbose) {
  #Simulate a single sample
  if (is.list(distrs)) distrsCur <- sample(distrs, 1)[[1]] else distrsCur <- distrs
  nSimReads <- N[which(N[,i] != 0),i]
  islandid <- names(nSimReads)
  curpis <- pis[,i]; names(curpis) <- rownames(pis)
  sim.r <- simReads(writeBam=0, seed=seed, islandid=islandid, nSimReads=nSimReads, pis=curpis, rl=as.integer(readLength), distrs=distrsCur, genomeDB=genomeDB, verbose=FALSE)
  zero <- which(N[,i]==0)
  v.list <- vector("list", length(zero))
  names(v.list) <-  rownames(N)[zero]
  sim.r@counts[[1]] <- c(v.list, sim.r@counts[[1]])
  exp.sim <- calcExp(distrs=distrsCur, genomeDB=genomeDB, pc=sim.r, readLength=readLength)
  ans <- list(explCnts=fData(exp.sim)[featnames,'explCnts',drop=FALSE], exp=exprs(exp.sim)[featnames,,drop=FALSE])
  if (verbose) cat('.')
  return(ans)
}


rowVar <- function (x, ...) { return((rowMeans(x^2, ...) - rowMeans(x, ...)^2) * ncol(x)/(ncol(x) - 1)) }



simnewsamplesNoisyObs <- function(fit,groupsnew,sel,x,groups,sigma2ErrorObs) {
# Simulate parameters from the posterior of NNGV under the presence of noisy observations, and new (non-noisy) observations from the posterior predictive
# Model
#  xhat_ij ~ N(x_ij, sigma2Error_i)
#  x_ij ~ N(mu_ik, sigma2Indiv_i), where k is group of observation j
#  mu_ik ~ Normal(mu0, sigma02)
#  phi_i = sigma2Error + sigma2Indiv ~ IGamma
#  sigma2Error ~ Gamma(a, l) where (a,l) are known (set via Method of Moments from sigma2ErrorObs)
#
# Equivalent marginal model for xhat_ij
#  xhat_ij ~ N(mu_ik, phi_i) where phi_i = sigma2Error_i + sigma2Indiv_i
#
# Returns posterior draws for (mu, phi, sigma2Indiv) and posterior predictive for x (the non-noisy version x rather than the noisy xhat)
if (is(x, "exprSet") | is(x,"ExpressionSet")) {
  if (is.character(groups) && length(groups)==1) { groups <- as.factor(pData(x)[, groups]) }
  x <- exprs(x)
} else if (!is(x,"data.frame") & !is(x,"matrix")) { stop("x must be an exprSet, data.frame or matrix") } 
patterns <- fit$patterns
groupsr <- gaga:::groups2int(groups,patterns); K <- as.integer(max(groupsr)+1)
groupsnewr <- gaga:::groups2int(groupsnew,patterns)
v <- fit$pp

# Checks
if ((max(groupsnewr)>max(groupsr)) | (min(groupsnewr)<min(groupsr))) stop('Groups indicated in groupsnew do not match with those indicated in groups')
if (ncol(x)!=length(groups)) stop('length(groups) must be equal to the number of columns in x')
if (!is.matrix(v)) stop('Argument v must be a matrix')
if (ncol(v)!=nrow(patterns)) stop('Argument v must have as many columns as rows has patterns')
if (nrow(x)!=nrow(v)) stop('Arguments x and v must have the same number of rows')
if (ncol(patterns)!=K) stop('patterns must have number of columns equal to the number of distinct elements in groups')
for (i in 1:nrow(patterns)) { patterns[i,] <- as.integer(as.integer(as.factor(patterns[i,]))-1) }
par <- getpar(fit)

sigmanew <- double(nrow(x))
munew <- matrix(NA,nrow=nrow(x),ncol=length(unique(groupsnewr))); colnames(munew) <- colnames(patterns)
xnew <- matrix(NA,nrow=nrow(x),ncol=length(groupsnewr))

#Draw delta
if (nrow(patterns)>1) {
  for (j in 2:ncol(v)) v[,j] <- v[,j]+v[,j-1]
  u <- runif(nrow(x)); d <- apply(u<v,1,function(z) which(z)[1])
} else {
  d <- rep(1,nrow(x))
}

#Draw sigma
a.sigma <- .5*(par['v0']+ncol(x))
b.sigma <- rep(.5*par['v0']*par['sigma02'], nrow(x))
simpat <- unique(d)
for (k in 1:length(simpat)) {
  rowsel <- d==k
  curgroups <- unique(patterns[k,])
  for (j in 1:length(curgroups)) {
    groupids <- colnames(patterns)[patterns[k,]==curgroups[j]]
    colsel <- groups %in% groupids
    xbar <- rowMeans(x[rowsel,colsel])
    b.sigma[rowsel] <- b.sigma[rowsel] + .5*rowSums((x[rowsel,colsel]-xbar)^2)
  }
}
sigmanew <- 1/sqrt(rgamma(nrow(x),a.sigma,b.sigma))

#Draw mu
for (k in 1:length(simpat)) {
  rowsel <- d==k
  curgroups <- unique(patterns[k,])
  for (j in 1:length(curgroups)) {
    groupids <- colnames(patterns)[patterns[k,]==curgroups[j]]
    colsel <- groups %in% groupids; ncolsel <- sum(colsel)
    xbar <- rowMeans(x[rowsel,colsel])
    v <- 1/(ncolsel/sigmanew[rowsel]^2 + 1/par['tau02'])
    m <- (ncolsel*xbar/sigmanew[rowsel]^2 + par['mu0']/par['tau02']) * v
    munew[rowsel, groupids] <- rnorm(sum(rowsel),m,sd=sqrt(v))
  }
}


#Draw x
sigmaIndiv <- sqrt(simEstError(sigma2=sigmanew^2, sigma2ErrorObs=sigma2ErrorObs)$sigma2Indiv)
design <- as.matrix(model.matrix(~ -1 + factor(groupsnewr)))
xnew <- t(design %*% t(munew)) + matrix(rnorm(nrow(x)*length(groupsnewr),0,sd=sigmaIndiv),nrow=nrow(xnew))

#Create ExpressionSet
metadata <- data.frame(labelDescription='Group that each array belongs to',row.names='group')
pheno <- new("AnnotatedDataFrame", data=data.frame(group=groupsnew), dimLabels=c("rowNames", "columnNames"), varMetadata=metadata)
sampleNames(pheno) <- paste("Array",1:nrow(pheno))
#
metadata <- data.frame(labelDescription=c('Expression pattern','SD for noisy obs','SD for non-noisy obs',paste('mean',colnames(patterns))),row.names=c('d','sigma','sigmaIndiv',paste('mean',1:ncol(patterns),sep='.')))
fdata <- data.frame(d,sigmanew,sigmaIndiv,munew)
fdata <- new("AnnotatedDataFrame", data=fdata,varMetadata=metadata)
sampleNames(fdata) <- paste("Gene",1:nrow(fdata))
varLabels(fdata) <- c('d','sigma','sigmaIndiv',paste('mean',1:ncol(patterns),sep='.'))
#
experimentData <- new("MIAME", title = "Dataset simulated with simnewsamples", abstract = "Expression data simulated from the posterior predictive distribution of a normal-normal model with generalized variances (see emfit in package EBarrays), under the assumption that input observations are noisy rather than exact. Posterior predictive observations are for non-noisy observations.")
#ans <- data.frame(xnew)
colnames(xnew) <- paste('Array',1:nrow(pheno)); rownames(xnew) <- paste('Gene',1:nrow(fdata))
xnew <- new("ExpressionSet", phenoData=pheno, featureData=fdata, exprs = xnew, experimentData = experimentData)
return(xnew)
}



simEstError <- function(sigma2, sigma2ErrorObs) {
  #Draw from posterior of simulation error variance sigma2Error ~ Gamma(a,l) I(sigma2Error <= sigma2)
  # Input
  # - sigma2: total variance is sigma2 = sigma2Error + sigma2Indiv, where sigma2Indiv is variation between individuals
  # - sigma2ErrorObs: posterior variance of estimates in pilot sample, used to set (a,l) by method of moments
  m <- rowMeans(sigma2ErrorObs); v <- rowVar(sigma2ErrorObs)
  l <- m/v; a <- m^2/v
  #logp <- pgamma(sigma2, a, l, log.p=TRUE)
  #logu <- -rtruncexp(length(logp), -logp)
  #logu[logu < -1e10] <- -1e10
  #sigma2Error <- qgamma(logu, a, l, log.p=TRUE)
  p <- pgamma(sigma2, a, l)
  u <- runif(length(sigma2), 0, p)
  sigma2Error <- qgamma(u, a, l)
  sel <- u==0; sigma2Error[sel] <- runif(sum(sel), 0, sigma2[sel]) #if round-off errors present, draw from uniform
  sel <- sigma2Error>sigma2; sigma2Error[sel] <- sigma2[sel]  #prevent round-off errors
  ans <- data.frame(sigma2Error=sigma2Error, sigma2Indiv=sigma2-sigma2Error)
  return(ans)
}


rtruncexp <- function(n, lower.trunc) {
  #Simulate n draws from lower-truncated exponential
  p <- pexp(lower.trunc)
  u <- runif(n, p, 1)
  qexp(u)
}
