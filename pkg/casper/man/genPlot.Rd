\name{genPlot}
\alias{genPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Plot exon structure and aligned reads for a given gene
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Plot exon structure for all transcripts in a given gene and aligned reads
}
\usage{
genPlot(goi, genomeDB, reads, exp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{goi}{
    %%     ~~Describe \code{goi} here~~
    ENTREZ id of gene of interest
}
  \item{genomeDB}{
    %%     ~~Describe \code{genomeDB} here~~
    List of annotations produced with the "procGenome" function
}
  \item{reads}{
    %%     ~~Describe \code{reads} here~~
    RangedData object of aligned reads or fragments 
}
  \item{exp}{
    %%     ~~Describe \code{exp} here~~
    ExpressionSet object with expression values
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{gene}{IRangesList object with one IRanges per transcript and
    its exons}
  \item{exp}{Named integer vector with transcript expression for the
    gene of interest}
%% ...

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
