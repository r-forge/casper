

modelPrior <- function(genomeDB, maxExons=40, smooth=TRUE, verbose=TRUE) {
  if (class(genomeDB) != "knownGenome") stop("genomeDB must be of class 'knownGenome'")

  if (verbose) cat("Counting number of annotated transcripts per gene... ")
  #CAMILLE: TO DO
  # - Compute table txsPerGene, which counts the nb of annotated transcripts per gene with 1,2... exons.
  #   rows correspond to nb of exons in the gene, columns to nb of annotated transcripts
  if (verbose) cat("Done.\n")

  if (verbose) cat("Counting number of exons contained in each variant... ")  
  #CAMILLE: TO DO
  # - Compute table exonsPerGene, which counts the nb of exons contained in a variant for genes with 1,2,... exons.
  #   rows correspond to nb of exons in the gene, columns to nb of exons
  if (verbose) cat("Done.\n")

  if (verbose) cat("Estimating parameters... ")
  nvarPrior <- nbVariantsDistrib(txsPerGene,maxExons=maxExons) #zero-truncated Negative Binomial fit
  nexonPrior <- nbExonsDistrib(exonsPerGene,maxExons=maxExons) #Beta-Binomial fit
  if (verbose) cat("Done.\n")
  return(list(nvarPrior=nvarPrior,nexonPrior=nexonPrior))
}

######################################################################
## PRIOR DISTRIBUTION FOR NUMBER OF EXONS
######################################################################

dbetabin <- function(x,n,alpha,beta) exp(lchoose(n,x) + lbeta(x+alpha,n-x+beta) - lbeta(alpha,beta))

nbExonsDistrib <- function(tab,maxExons=40,smooth=TRUE) {
  #Model nb exons per variant E_d as (E_d-1) ~ Beta-Binomial(E,alpha,beta), where E=nb exons in the gene.
  # Input
  # - tab: table counting nb exons (rows) and nb annotated transcripts (columns). rownames and colnames must be specified.
  # - maxExons: counts for genes with >maxExons are merged.
  # - smooth: if set to TRUE estimates for genes with 10 up to maxExons are smoothed via cross-validated gam. Improves stability.
  # Ouput
  # - bbpar: estimated Beta-Binomial parameters. For genes with > maxExons, the estimate for maxExons is carried over
  # - obs: list. Element i is the empirical distribution of nb variants for genes with i exons
  # - pred: list with predicted distributions

  require(VGAM)

  sel <- as.numeric(rownames(tab))>maxExons
  extrapolate <- tab[sel,]
  tab <- tab[!sel,]

  #Estimate parameters
  bbpar <- matrix(NA,nrow=nrow(tab)+nrow(extrapolate),ncol=2); colnames(bbpar) <- c('alpha','beta')
  obs <- pred <- vector("list",nrow(tab)+nrow(extrapolate))
  rownames(bbpar) <- names(obs) <- names(pred) <- c(rownames(tab),rownames(extrapolate))
  bbpar[1,] <- c(0,0)
  for (i in nrow(tab):2) {
    y <- tab[i,tab[i,]>0]
    if (i==2) {
      bbpar[i,1] <- sum((as.numeric(names(y))-1) * y / sum(y)) * sum(bbpar[i+1,])
      bbpar[i,2] <- sum(bbpar[i+1,]) - bbpar[i,1]
    } else {
      ydf <- rep(as.numeric(names(y))-1,y) #start at 0
      ydf <- data.frame(succ=ydf,fail=i-ydf-1)
      fit <- vglm(cbind(succ, fail) ~ 1, betabinomial.ab, data=ydf, trace=FALSE)
      #fit <- vglm(cbind(succ,fail) ~ 1, family=betabinomial, data=y, trace=FALSE)
      bbpar[i,] <- Coef(fit)
    }
  }

  #smooth parameter estimates for genes with >=10 exons
  if (smooth==TRUE) {
    require(mgcv)
    m <- bbpar[2:nrow(tab),1]/rowSums(bbpar)[2:nrow(tab)]
    m <- data.frame(logitm= log(m/(1-m)), E= as.numeric(rownames(bbpar)[2:nrow(tab)]))
    fit <- gam(logitm ~ s(E, sp= -1), data=m)
    msmooth <- 1/(1+exp(-predict(fit)))
    b <- data.frame(b=bbpar[2:nrow(tab),2], E= as.numeric(rownames(bbpar)[2:nrow(tab)]))
    fit <- gam(b ~ s(E, sp= -1), data=b, family=gaussian(link="log"))
    bsmooth <- exp(predict(fit))
    asmooth <- bsmooth*msmooth/(1-msmooth)
     
    sel <- as.character(m$E[m$E>=10])
    bbpar[sel,1] <- asmooth[sel]
    bbpar[sel,2] <- bsmooth[sel]
  }

  bbpar[(nrow(tab)+1):nrow(bbpar),] <- rep(bbpar[nrow(tab),],each=nrow(bbpar)-nrow(tab))
  
  #predicted frequencies
  pred[[1]] <- obs[[1]] <- tab['1',]
  names(pred[[1]]) <- names(obs[[1]]) <- '1'
  nkeep <- nrow(tab)+1
  tab <- rbind(tab,extrapolate)
  for (n in names(pred)[-1]) {
    i <- as.numeric(n)
    y <- tab[n,tab[n,]>0]
    pred[[n]] <- sum(y)*dbetabin(0:(i-1),n=i-1,alpha=bbpar[n,1],beta=bbpar[n,2])
    names(pred[[n]]) <- 1:i
    obs[[n]] <- y[names(pred[[n]])]
    obs[[n]][is.na(obs[[n]])] <- 0
    names(obs[[n]]) <- names(pred[[n]])
  }

  return(list(bbpar=bbpar[1:nkeep,],obs=obs,pred=pred))
}


######################################################################
## PRIOR DISTRIBUTION FOR NUMBER OF VARIANTS
######################################################################

nbVariantsDistrib <- function(tab,maxExons=40) {
  # Input
  # - tab: table counting nb exons (rows) and nb annotated transcripts (columns). rownames and colnames must be specified.
  # - maxExons: counts for genes with >maxExons are merged.
  # Ouput
  # - nbpar: estimated Negative Binomial parameters (size: nb coin flips; prob: success prob)
  # - obs: list. Element i is the empirical distribution of nb variants for genes with i exons
  # - pred: list with predicted distributions
  
  sel <- as.numeric(rownames(tab))>maxExons
  tab <- rbind(tab[!sel,],colSums(tab[sel,,drop=FALSE]))
  rownames(tab)[maxExons+1] <- as.character(maxExons+1)
  
  n <- as.numeric(colnames(tab))
  maxtxs <- log2(2^(as.numeric(rownames(tab))) - 1) #nb variants cannot be > 2^p -1
  nbpar <- matrix(NA,nrow=nrow(tab),ncol=2); colnames(nbpar) <- c('prob','size')
  obs <- pred <- vector("list",nrow(tab))
  names(maxtxs) <- rownames(nbpar) <- names(obs) <- names(pred) <- rownames(tab)
  nbpar[1,] <- c(0,0)
  ichar <- rownames(tab)
  for (i in 2:nrow(tab)) {
    y <- tab[ichar[i],log2(n)<=maxtxs[ichar[i]]]
    if (length(y)>1) {
      fit <- getNBinomParams(y,components=1)
      nbpar[ichar[i],] <- fit[c('p','s1')]
    } else {
      nbpar[ichar[i],] <- nbpar[ichar[i-1],]
    }
    ymax <- max(as.numeric(names(y)))
    p <- dnbinom(1:ymax,size=nbpar[ichar[i],'size'],prob=nbpar[ichar[i],'prob'])/dnbinom(0,size=nbpar[ichar[i],'size'],prob=nbpar[ichar[i],'prob'])
    pred[[ichar[i]]] <- sum(y)*p/sum(p)
    names(pred[[ichar[i]]]) <- 1:ymax
    obs[[ichar[i]]] <- y[names(pred[[ichar[i]]])]
    obs[[ichar[i]]][is.na(obs[[ichar[i]]])] <- 0
    names(obs[[ichar[i]]]) <- names(pred[[ichar[i]]])
  }
  return(list(nbpar=nbpar,obs=obs,pred=pred))
}


#Estimate parameters of zero-truncated Negative Binomial (set components>1 to fit mixture)
#Input x is a vector with counts. names(x) must be specified.
# e.g. x= c(100,50,10,1); names(x)= c('1','2','3','5')
getNBinomParams <- function(x,mc.cores=1,components=1) {
  logit <- function(x) { return(log(x/(1-x))) }
  ilogit <- function(x) { return(1/(1+exp(-x))) }
  myLikelihood <- function(params,x,components) {
    if (components==1) {
      size1 <- params[1]; prob1 <- params[2]
      f1 <- dnbinom(as.numeric(names(x)),size1,prob1)
      f1.sum <- sum(f1)
      if (is.nan(f1.sum) | f1.sum==0) {
        ans <- Inf
      } else {
        q <- f1/f1.sum
        ans <-  -dmultinom(x=x,prob=q,log=TRUE)
      }
    } else if (components==2) {
      size1 <- params[1]; size2 <- params[2]; prob1 <- params[3]; prob2 <- params[4]; w <- params[5]
      f1 <- dnbinom(as.numeric(names(x)),size1,prob1)
      f2 <- dnbinom(as.numeric(names(x)),size2,prob2)
      f1.sum <- sum(f1,na.rm=TRUE)
      f2.sum <- sum(f2,na.rm=TRUE)
      if (is.nan(w) | is.nan(f1.sum) | is.nan(f2.sum) | is.na(w) | f1.sum==0 | f2.sum==0) {
        ans <- Inf
      } else if (w<0.5) {
        ans <- Inf
      } else {
        q <- (w * f1/f1.sum) + ((1-w) * f2/f2.sum)
        ans <-  -dmultinom(x=x,prob=q,log=TRUE)
      }
    } else if (components==3) {
      size1 <- params[1]; size2 <- params[2]; size3 <- params[3]; prob1 <- params[4]; prob2 <- params[5]; prob3 <- params[6]; w1 <- params[7]; w2 <- params[8]
      f1 <- dnbinom(as.numeric(names(x)),size1,prob1)
      f2 <- dnbinom(as.numeric(names(x)),size2,prob2)
      f3 <- dnbinom(as.numeric(names(x)),size3,prob3)      
      f1.sum <- sum(f1,na.rm=TRUE)
      f2.sum <- sum(f2,na.rm=TRUE)
      f3.sum <- sum(f3,na.rm=TRUE)
      if (is.nan(w1) | is.nan(w2) | is.nan(f1.sum) | is.nan(f2.sum) | is.nan(f3.sum) | is.na(w1) | is.na(w2) | f1.sum==0 | f2.sum==0 | f3.sum==0) {
        ans <- Inf
      } else if (w1<0.5 | (w1+w2>=1)) {
        ans <- Inf
      } else {
        q <- (w1 * f1/f1.sum) + (w2 * f2/f2.sum) + ((1 - w1 - w2) * f3/f3.sum)
        ans <-  -dmultinom(x=x,prob=q,log=TRUE)
      }
    } else if (components==4) {
      size1 <- params[1]; size2 <- params[2]; size3 <- params[3]; size4 <- params[4]
      prob1 <- params[5]; prob2 <- params[6]; prob3 <- params[7]; prob4 <- params[8]
      w1 <- params[9]; w2 <- params[10]; w3 <- params[11]
      f1 <- dnbinom(as.numeric(names(x)),size1,prob1)
      f2 <- dnbinom(as.numeric(names(x)),size2,prob2)
      f3 <- dnbinom(as.numeric(names(x)),size3,prob3)
      f4 <- dnbinom(as.numeric(names(x)),size4,prob4)              
      f1.sum <- sum(f1,na.rm=TRUE)
      f2.sum <- sum(f2,na.rm=TRUE)
      f3.sum <- sum(f3,na.rm=TRUE)
      f4.sum <- sum(f4,na.rm=TRUE)
      if (is.nan(w1) | is.nan(w2) | is.nan(w3) | is.nan(f1.sum) | is.nan(f2.sum) | is.nan(f3.sum) | is.nan(f4.sum) | is.na(w1) | is.na(w2) | is.na(w3) | f1.sum==0 | f2.sum==0 | f3.sum==0 | f4.sum==0) {
        ans <- Inf          
      } else if (w1<0.5 | (w1+w2+w3>=1)) {
        ans <- Inf
      } else {
        q <- (w1 * f1/f1.sum) + (w2 * f2/f2.sum) + (w3 * f3/f3.sum) + ((1 - w1 - w2 - w3) * f4/f4.sum)
        ans <-  -dmultinom(x=x,prob=q,log=TRUE)
      }
    }
    ans
  }
  fopt <- function(params,x,giveMe,components) {
    #reverse log and logit
    params[grepl('[pw]',names(params))] <- ilogit(params[grepl('[pw]',names(params))])
    params[grepl('s',names(params))] <- exp(params[grepl('s',names(params))])
    #run mylikelihood
    ans <- myLikelihood(params,x,components)
    ans
  }
  mynlminb <- function(params,numexp,numtrial,giveMe='objective',components) {
    #add size as parameter
    if (components==1) {
      s1 <- numtrial * params[1] / numexp
      names(params) <- c('p')
      params <- c(s1=s1,params)
    } else if (components==2) {
      s1 <- numtrial * params[1] / numexp; s2 <- numtrial * params[2] / numexp
      names(params) <- c('p1','p2','w')        
      params <- c(s1=s1,s2=s2,params)
    } else if (components==3) {
      s1 <- numtrial * params[1] / numexp; s2 <- numtrial * params[2] / numexp; s3 <- numtrial * params[3] / numexp
      names(params) <- c('p1','p2','p3','w1','w2') 
      params <- c(s1=s1,s2=s2,s3=s3,params)
    } else if (components==4) {
      s1 <- numtrial * params[1] / numexp; s2 <- numtrial * params[2] / numexp; s3 <- numtrial * params[3] / numexp; s4 <- numtrial * params[4] / numexp
      names(params) <- c('p1','p2','p3','p4','w1','w2','w3') 
      params <- c(s1=s1,s2=s2,s3=s3,s4=s4,params)
    }
    #log and logit
    params[grepl('[pw]',names(params))] <- logit(params[grepl('[pw]',names(params))])
    params[grepl('s',names(params))] <- log(params[grepl('s',names(params))])
    #run optimizer
    ans <- nlminb(start=params,function(y) fopt(y,x,giveMe,components))
    ans[[giveMe]]
  }
  df2list <- function(y) unclass(as.data.frame(t(unique(do.call(rbind,lapply(p,function(x) y))))))
  p <- seq(0.1,0.9,length=4); w1 <- seq(0.6,1,length=3); w2 <- seq(0.1,0.4,length=3)
  tmp <- do.call(rbind,lapply(p,function(x) cbind(rep(x,length(p)),p)))
  tmp <- do.call(rbind,lapply(p,function(x) cbind(tmp,rep(x,length(p)))))
  tmp <- do.call(rbind,lapply(p,function(x) cbind(tmp,rep(x,length(p)))))    
  tmp <- do.call(rbind,lapply(w1,function(x) cbind(tmp,rep(x,nrow(tmp)))))
  tmp <- do.call(rbind,lapply(w2,function(x) cbind(tmp,rep(x,nrow(tmp)))))
  tmp <- do.call(rbind,lapply(w2,function(x) cbind(tmp,rep(x,nrow(tmp)))))
  tmp <- tmp[tmp[,1]!=tmp[,2] & tmp[,1]!=tmp[,3] & tmp[,1]!=tmp[,4] & tmp[,2]!=tmp[,3] & tmp[,2]!=tmp[,4] & tmp[,3]!=tmp[,4] & (tmp[,4]+tmp[,5])<=1 & (tmp[,5]+tmp[,6])<=1 & (tmp[,5]+tmp[,6]+tmp[,7])<=1,]
  numexp <- sum(x)
  numtrial <- sum(x * as.numeric(names(x)))
  if (components==0) {
    tmp.1.list <- unclass(as.data.frame(t(unique(tmp[,1,drop=FALSE]))))
    tmp.2.list <- unclass(as.data.frame(t(unique(tmp[,c(1,2,5)]))))
    if (length(use)>8) tmp.3.list <- unclass(as.data.frame(t(unique(tmp[,c(1,2,3,5,6)]))))
    if (length(use)>11) tmp.4.list <- unclass(as.data.frame(t(unique(tmp))))
    if (mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        #
        myfun <- function(idx) {
          ans <- vector('list',length(idx))
          for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.1.list[[idx[i]]],numexp,numtrial,'objective',components=1)
          return(ans)
        }
        objec.1 <- unlist(multicore::pvec(1:length(tmp.1.list),myfun,mc.cores=ifelse(length(tmp.1.list)<=mc.cores,round(length(tmp.1.list)/2),mc.cores)))
        params.1 <- mynlminb(tmp.1.list[[which(objec.1==min(objec.1))[1]]],numexp,numtrial,'par',components=1)
        toadd <- df2list(cbind(rep(ilogit(params.1[2]),length(p)),p,1))
        tmp.2.list[(length(tmp.2.list)+1):(length(tmp.2.list)+length(toadd))] <- toadd
        #
        myfun <- function(idx) {
          ans <- vector('list',length(idx))
          for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.2.list[[idx[i]]],numexp,numtrial,'objective',components=2)
          return(ans)
        }
        objec.2 <- unlist(multicore::pvec(1:length(tmp.2.list),myfun,mc.cores=ifelse(length(tmp.2.list)<=mc.cores,round(length(tmp.2.list)/2),mc.cores)))
        params.2 <- mynlminb(tmp.2.list[[which(objec.2==min(objec.2))[1]]],numexp,numtrial,'par',components=2)
        toadd <- df2list(cbind(rep(ilogit(params.2[3]),length(p)),rep(ilogit(params.2[4]),length(p)),p,rep(ilogit(params.2[5]),length(p)),1-ilogit(params.2[5])))
        if (length(use)>8) {
          tmp.3.list[(length(tmp.3.list)+1):(length(tmp.3.list)+length(toadd))] <- toadd
          #
          myfun <- function(idx) {
            ans <- vector('list',length(idx))
            for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.3.list[[idx[i]]],numexp,numtrial,'objective',components=3)
            return(ans)
          }
          objec.3 <- unlist(multicore::pvec(1:length(tmp.3.list),myfun,mc.cores=ifelse(length(tmp.3.list)<=mc.cores,round(length(tmp.3.list)/2),mc.cores)))
          params.3 <- mynlminb(tmp.3.list[[which(objec.3==min(objec.3))[1]]],numexp,numtrial,'par',components=3)
          toadd <- df2list(cbind(rep(ilogit(params.3[4]),length(p)),rep(ilogit(params.3[5]),length(p)),rep(ilogit(params.3[6]),length(p)),p,rep(ilogit(params.3[7]),length(p)),rep(ilogit(params.3[8]),length(p)),1-ilogit(params.3[7])-ilogit(params.3[8])))
        }
        if (length(use)>11) {          
          tmp.4.list[(length(tmp.4.list)+1):(length(tmp.4.list)+length(toadd))] <- toadd
          myfun <- function(idx) {
            ans <- vector('list',length(idx))
            for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.4.list[[idx[i]]],numexp,numtrial,'objective',components=4)
            return(ans)
          }
          objec.4 <- unlist(multicore::pvec(1:length(tmp.4.list),myfun,mc.cores=ifelse(length(tmp.4.list)<=mc.cores,round(length(tmp.4.list)/2),mc.cores)))
          params.4 <- mynlminb(tmp.4.list[[which(objec.4==min(objec.4))[1]]],numexp,numtrial,'par',components=4)
        }
      } else stop('multicore library has not been loaded!')
    } else {
      objec.1 <- unlist(lapply(tmp.1.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components=1)))
      params.1 <- mynlminb(tmp.1.list[[which(objec.1==min(objec.1))[1]]],numexp,numtrial,'par',components=1)
      toadd <- df2list(cbind(rep(ilogit(params.1[2]),length(p)),p,1))
      tmp.2.list[(length(tmp.2.list)+1):(length(tmp.2.list)+length(toadd))] <- toadd
      objec.2 <- unlist(lapply(tmp.2.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components=2)))
      params.2 <- mynlminb(tmp.2.list[[which(objec.2==min(objec.2))[1]]],numexp,numtrial,'par',components=2)
      toadd <- df2list(cbind(rep(ilogit(params.2[3]),length(p)),rep(ilogit(params.2[4]),length(p)),p,rep(ilogit(params.2[5]),length(p)),1-ilogit(params.2[5])))
      if (length(use)>8) {
        tmp.3.list[(length(tmp.3.list)+1):(length(tmp.3.list)+length(toadd))] <- toadd
        objec.3 <- unlist(lapply(tmp.3.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components=3)))
        params.3 <- mynlminb(tmp.3.list[[which(objec.3==min(objec.3))[1]]],numexp,numtrial,'par',components=3)
        toadd <- df2list(cbind(rep(ilogit(params.3[4]),length(p)),rep(ilogit(params.3[5]),length(p)),rep(ilogit(params.3[6]),length(p)),p,rep(ilogit(params.3[7]),length(p)),rep(ilogit(params.3[8]),length(p)),1-ilogit(params.3[7])-ilogit(params.3[8])))
      }
      if (length(use)>11) {
        tmp.4.list[(length(tmp.4.list)+1):(length(tmp.4.list)+length(toadd))] <- toadd
        objec.4 <- unlist(lapply(tmp.4.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components=4)))
        params.4 <- mynlminb(tmp.4.list[[which(objec.4==min(objec.4))[1]]],numexp,numtrial,'par',components=4)
      }
    }
    #Bayesian information criterion.
    bic.1 <- -2 * (-min(objec.1)) + 2 * log(sum(x)) #objec contains the -log(likelihood) 
    bic.2 <- -2 * (-min(objec.2)) + 5 * log(sum(x))
    if (length(use)>8) bic.3 <- -2 * (-min(objec.3)) + 8 * log(sum(x))
    if (length(use)>11) bic.4 <- -2 * (-min(objec.4)) + 11 * log(sum(x))
    if (length(use)>11) {
      bic <- c(bic.1,bic.2,bic.3,bic.4)
    } else if (length(use)>8) {
      bic <- c(bic.1,bic.2,bic.3)        
    } else {
      bic <- c(bic.1,bic.2)
    }
    components <- which(bic==min(bic))
    if (length(use)>11) {
      ans <- params.list <- list(params.1,params.2,params.3,params.4)[[components]]
    } else if (length(use)>8) {
      ans <- params.list <- list(params.1,params.2,params.3)[[components]]
    } else {
      ans <- params.list <- list(params.1,params.2)[[components]]
    }
  } else {
    if (components==1) {
      tmp <- unique(tmp[,1,drop=FALSE])
    } else if (components==2) {
      tmp <- unique(tmp[,c(1,2,5)])
    } else if (components==3) {
      tmp <- unique(tmp[,c(1,2,3,5,6)])
    }
    tmp.list <- unclass(as.data.frame(t(tmp)))
    if (mc.cores>1) {
      if ('multicore' %in% loadedNamespaces()) {
        myfun <- function(idx) {
          ans <- vector('list',length(idx))
          for (i in 1:length(idx)) ans[[i]] <- mynlminb(params=tmp.list[[idx[i]]],numexp,numtrial,'objective',components)
          return(ans)
        }
        objec <- unlist(multicore::pvec(1:length(tmp.list),myfun,mc.cores=ifelse(length(tmp.list)<=mc.cores,round(length(tmp.list)/2),mc.cores)))
      } else stop('multicore library has not been loaded!')
    } else {
      objec <- unlist(lapply(tmp.list,function(x) mynlminb(params=x,numexp,numtrial,'objective',components)))
    }
    ans <- mynlminb(tmp.list[[which(objec==min(objec))[1]]],numexp,numtrial,'par',components)
  }
  ans[grepl('[pw]',names(ans))] <- ilogit(ans[grepl('[pw]',names(ans))])
  ans[grepl('s',names(ans))] <- exp(ans[grepl('s',names(ans))])
  ans
}
