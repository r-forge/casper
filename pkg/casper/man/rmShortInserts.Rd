\name{rmShortInserts}
\alias{rmShortInserts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Remove reads with short insert sizes from imported BAM files. 
}
\description{
  In paired-end experiments short inserts indicate that the two ends are very close to each other,
  which can indicate RNA degradation.
  The insert size is stored in objects imported with \code{scanBam} in the element named \code{isize}.
}
\usage{
rmShortInserts(bam, isizeMin=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bam}{ Object with aligned reads, as returned by \code{scanBam}  }
  \item{isizeMin}{ Reads with insert size smaller than \code{isizeMin} will be removed. }
}
\value{
Named list, in the same format as that returned by \code{scanBam}.
}
\author{
David Rossell
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~paired-end sequencing }
\keyword{ ~SAM/BAM }% __ONLY ONE__ keyword per line
