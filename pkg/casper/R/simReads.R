simReads <- function(selIslands, nSimReads, pis, rl, seed, writeBam, distrs, genomeDB, chrlen=NULL, repSims=FALSE, bamFile=NULL, stranded=FALSE, samtoolsPath=NULL, mc.cores=1){
  ##### Simulate reads
  if(writeBam==1) {
    if(is.null(chrlen)){
      fn <- paste("http://hgdownload.cse.ucsc.edu/goldenPath/",genomeDB@genomeVersion,"/database/chromInfo.txt.gz", sep="")
      dfile="tmp.gz"
      cat(paste("Downloading chromosome lengths from UCSC for genome ", genomeDB@genomeVersion, " (file: ",fn,"\nUsing ", dfile, " as destination file (deleted afterwards)\n", sep=""))
      dn <- try(download.file(fn,destfile=dfile, quiet=T), silent=T)
      if(class(dn)=="try-error") stop("No chromInfo file found at htt://hgdownload.cse.ucsc.edu/goldenPath/ for required genome. Please provide chrlen vector with reference sequence chromosomes lengths")
      chrlen <- readLines(gzfile(dfile))
      chrlen <- strsplit(chrlen, split="\t")
      chrlen <- lapply(chrlen, '[', 1:2)
      chrlen <- do.call(rbind, chrlen)
      namchr <- chrlen[,1]
      chrlen <- as.numeric(chrlen[,2])
      names(chrlen) <- namchr
      system(paste("rm ", dfile))
    }
    if(!grepl(".bam$", bamFile)) stop("bamFile must be a valid bam file name (end with .bam)")
    if(is.null(samtoolsPath)) stop("samtoolsPath must be specified when writeBam=TRUE")
  }
  tmp <- sub(" ", ".", Sys.time())
  lr_file <- paste(bamFile, ".lr.", tmp, ".sam", sep="")
  rr_file <- paste(bamFile, ".rr.", tmp, ".sam", sep="")
  sims <- casperSim(genomeDB=genomeDB, distrs=distrs, nSimReads=nSimReads, pis=pis, selIslands=selIslands, lr_file=lr_file, rr_file=rr_file, rl=rl, chrlen, seed, bam=writeBam)
  if(writeBam==1){
    finalBamFile=paste(rr_file, ".fin.bam", sep="")
    finalBamSorted=sub(".bam", "", bamFile)
    cat(paste("generating bam ", bamFile,"\n"))
    system(paste("grep -v '@' ", lr_file, "  >> ", rr_file, " && ", samtoolsPath," view -Sb ", rr_file, "  > ", finalBamFile, " && ", samtoolsPath," sort ", finalBamFile, " ", finalBamSorted," && ", samtoolsPath, " index ", bamFile, sep=""))
    system(paste("rm ", lr_file, " && rm ", rr_file, " && rm ", finalBamFile, sep=""))
  }
  cat("Splitting counts\n")
  counts <- splitPaths(sims$pc, genomeDB, mc.cores=mc.cores, stranded=stranded, geneid=selIslands)
  pc <- new("pathCounts", counts=counts, denovo=FALSE, stranded=stranded)
  if(repSims==1) {
    sims$pc <- NULL
    sims$Nsim <- NULL
    ans <- list(sims=sims, pc=pc)
  }
  else ans <- pc
  ans
}

splitPaths <- function(paths, DB, mc.cores, stranded, geneid){
  paths <- paths[grepl("^\\..*\\.$", names(paths))]
  sel <- strsplit(names(paths), split='-|\\.')
  sel1 <- lapply(sel, "[", 2)
  sel1 <- unlist(sel1)
  islands <- DB@islands[geneid]
  nislEx <- lapply(islands, length)
  nislEx <- rep(names(islands), unlist(nislEx))
  islEx <- unlist(lapply(islands, names))
  names(islEx) <- nislEx
  isl <- match(sel1, islEx)
  isl <- names(islEx)[isl]
  splCounts <- split(paths, isl)
  splCounts <- lapply(splCounts, function(x) x[grepl("-", names(x))])
  if(DB@denovo){
    sel <- sapply(sel, "[", -1)
    tmp <- split(sel, isl)
    if(mc.cores>1) {
      require(multicore)
      tmp1 <- multicore:::mclapply(names(tmp), function(x){
        n <- sapply(tmp[[x]], length)
        nn <- unlist(tmp[[x]])
        names(nn) <- rep(names(splCounts[[x]]), n)
        nnn <- nn %in% names(DB@islands[[x]])
        nnnn <- tapply(nnn, names(nn), all)
        splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
      }, mc.cores=mc.cores)
    } else {
      tmp1 <- lapply(names(tmp), function(x){
        n <- sapply(tmp[[x]], length)
        nn <- unlist(tmp[[x]])
        names(nn) <- rep(names(splCounts[[x]]), n)
        nnn <- nn %in% names(DB@islands[[x]])
        nnnn <- tapply(nnn, names(nn), all)
        splCounts[[x]][names(splCounts[[x]]) %in% names(nnnn)[nnnn]]
      })
    }
    names(tmp1) <- names(tmp)
    splCounts <- tmp1
  }
  ans <- vector(length(islands), mode='list')
  names(ans) <- names(islands)
  ans[names(splCounts)] <- splCounts
  if(!stranded) ans <- list(ans)
  else {
    plus <- ans[names(ans) %in% names(DB@islandStrand)[DB@islandStrand=='+']]
    minus <- ans[names(ans) %in% names(DB@islandStrand)[DB@islandStrand=='-']]
    ans <- list(minus=minus, plus=plus)
  }
  ans
}

casperSim <- function(genomeDB, distrs, nSimReads, pis, selIslands, lr_file=NULL, rr_file=NULL, rl, chrlen, seed, bam){
  cat("Formatting input\n")
  sel <- selIslands
  txs <- genomeDB@transcripts[sel]
  sel1 <- unlist(lapply(txs, is.null))
  txs <- txs[!sel1]
  island_strand <- genomeDB@islandStrand[sel][!sel1]
  island_strand[island_strand=="+"] <- 1
  island_strand[island_strand=="-"] <- -1
  tmp <- names(island_strand)
  island_strand <- as.numeric(island_strand)
  names(island_strand) <- tmp
  variant_num <- unlist(lapply(txs, length))
  starts <- start(genomeDB@exonsNI)
  names(starts) <- genomeDB@exonsNI$id
  ends <- end(genomeDB@exonsNI)
  names(ends) <- names(starts)
  chroms <- genomeDB@exonsNI$space
  names(chroms) <- genomeDB@exonsNI$id
  exon_num <- unlist(lapply(txs, function(x) lapply(x, length)))
  exs <- unlist(txs)
  exon_st <- starts[as.character(exs)]
  exon_end <- ends[as.character(exs)]
  exs <- unname(exs)
  sel2 <- unlist(lapply(lapply(txs, '[[', 1), '[', 1))
  chroms <- as.character(chroms[as.character(sel2)])
  tmps <- rep(1:length(exon_num), exon_num)
  widths <- exon_end-exon_st+1
  variant_len <- tapply(widths, tmps, sum)
  txs <- unlist(lapply(txs,names)) 
  if(!all(txs %in% names(pis))) stop("Wrong pis vector, some transcripts missing")
  if(!all(names(nSimReads) %in% selIslands)) stop("Wrong nSimReads vector, some islands missing")
  if(!all(nSimReads>0)) stop("nSimReads with zero entries")
  ge <- nSimReads[selIslands]
  ge <- as.integer(rep(0:(length(ge)-1), ge))
  ve=pis[txs]
  gs = island_strand
  vn=variant_num
  vl=variant_len
  en=exon_num
  es=exon_st
  ee=exon_end
  ei=exs
  ngenes=length(sel)
  ldv <- sample(as.numeric(names(distrs@lenDis)), p=distrs@lenDis/sum(distrs@lenDis), size=sum(nSimReads), replace=T)
  ldd <- as.integer(1)
  th <- seq(0, 1, len=10000)
  std <- distrs@stDis(th)
  std[1] <- 0
  std[length(th)] <- 1
  sdv <- approxfun(std[!is.na(std)], th[!is.na(std)])(th)
  sdd <- std
  if(bam) {
    write.sam.header(chrlen[names(chrlen) %in% unique(chroms)], lr_file, max(vn))
    write.sam.header(chrlen[names(chrlen) %in% unique(chroms)], rr_file, max(vn))
  }
  cat("Simulating fragments\n")
# insideBam deprecated, only in case it is useful in simulations, return from C all information written to the bam file
  insideBam=integer(0)
  ans <- .Call("casperSimC", ge, ve, vn, vl, en, es, ee, ei, ldv, ldd, sdv, sdd, rl, length(ge), gs, lr_file, rr_file, chroms, as.integer(seed), as.integer(bam), as.integer(insideBam))
  ans <- ans[1:7]
  names(ans[[7]]) <- ans[[6]]
  ans[[6]] <- NULL
  names(ans) <- c("varl", "st", "len", "abst", "strand", "pc")
  ans
}

write.sam.header <- function(chrlen, file, nvars){
  chr.line <- function(chr, len, file, append) cat("@SQ\tSN:", chr, "\tLN:", len, "\n", sep="", file=file, append=append)
  chr.line(names(chrlen)[1], chrlen[1], file, append=F)
  if(length(chrlen)>1) for(i in 2:length(chrlen)) chr.line(names(chrlen)[i], chrlen[i], file, append=T)
  cat("@PG\tID:simulation\tPN:simulation\tVN:0.0.1\n", file=file, append=T)
  for(i in 1:nvars) cat("@RG\tID:", i, "\tSM:", i, "\n", file=file, append=T)
}

