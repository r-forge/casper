\name{pathCounts}
\alias{pathCounts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Compute exon path counts
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Compute counts for exon paths visited by aligned reads
}
\usage{
pathCounts(reads, exons)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads}{
    %%     ~~Describe \code{reads} here~~
    RangeData object with aligned reads
}
  \item{exons}{
    %%     ~~Describe \code{exons} here~~
    RangedData object with non-overlapping exons
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  Named integer vector with counts of exon paths. Names are character
  strings built as ".exon1.exon2-exon3.exon4.", with dashes making the
  split between exons visited by left and right-end reads
  correspondingly. 
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
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
