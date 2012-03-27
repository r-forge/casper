\name{getDistrs}
\alias{getDistrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Compute fragment start and fragment length distributions
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Compute fragment start distributions by using reads
  aligned to genes with only one annotated variant.
  Estimate fragment length
  distribution using fragments aligned to long exons (>1000nt).
  Fragment length is defined as the distance between the start of the
  left-end read and the end of the right-end read.
}
\usage{
getDistrs(txs, exons, frags, nreads=4*10^6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txs}{
    GRanges object with known transcripts
}
  \item{exons}{
    GRangesList with annotated exons per transcript
}
  \item{frags}{
    List with reads from imported BAM file, as returned by
    \code{scanBam} or \code{rmShortInserts}.
  }
  \item{nreads}{
    To speed up computations, only the first \code{nreads} are used to
    obtain the estimates. The default value of 4 milions usually gives
    highly precise estimates.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  \item{lenDis}{Table with number of fragments with a given length}
  \item{stDis}{Numeric vector with relative fragment start positions}
}
\author{
Camille Stephan-Otto Attolini, David Rossell
}
\examples{
##Process genome
##hg19DB <- procGenome('hg19',mc.cores=6)

##Read BAM file
##what <- c('qname','flag','rname','strand','pos','qwidth','mapq','cigar','mrnm','mpos','isize')
##param <- ScanBamParam(scanBamFlag(isPaired=TRUE,hasUnmappedMate=FALSE),what=what)
##bam0 <- scanBam(file='~/Desktop/ranidata/test_rani/mapped_reads/accepted_hits.bam',param=param)[[1]]

##Remove short inserts
##bam0 <- rmShortInserts(bam0, isizeMin=150)

##Estimate fragment size distribution
##distrs <- getDistrs(txs=hg19DB$txs,exons=hg19DB$exons,bam=bam0,nreads=length(bam0$qname))

##Plot fragment size distribution
##n <- as.numeric(names(distrs$lenDis))
##plot(n,distrs$lenDis,type='l',xlab='Fragment length',ylab='Counts')

##Plot read start distribution
##hist(distrs$stDis,xlab='Read start',main='')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stats }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
