\name{getDistrs}
\alias{getDistrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Compute fragment start and fragment length distributions
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  This function calculates fragment start distributions for reads
  aligned to genes with only one annotated variant and fragment length
  distribution for fragments aligned to long exons (>1000nt)
}
\usage{
getDistrs(txs, exons, frags, mc.cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{txs}{
    %%     ~~Describe \code{txs} here~~
    GRanges object with known transcripts
}
  \item{exons}{
    %%     ~~Describe \code{exons} here~~
    GRangesList with annotated exons per transcript
}
  \item{frags}{
    %%     ~~Describe \code{frags} here~~
    RangedData object defined by start and end of whole fragments
}
  \item{mc.cores}{
    %%     ~~Describe \code{mc.cores} here~~
    Number of cores to use in parallel processing (multicore package
    required)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  \item{stDis}{Numeric vector with relative fragment start positions}
  \item{lenDis}{Table with fragment counts for each existing length}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Camille Stephan-Otto Attolini
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line