\name{procGenome}
\alias{procGenome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Download and format annotations for a given genome
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Download USCS annotations for a given genome
}
\usage{
procGenome(genome, mc.cores = mc.cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome}{
    %%     ~~Describe \code{genome} here~~
    Genome version from UCSC (e.g. "hg19", "dm3")
}
  \item{mc.cores}{
    %%     ~~Describe \code{mc.cores} here~~
    Number of cores to use in parallel processing (multicore package required)
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  This function generates all necessary annotation objects for subsequent
  functions.
}
\value{
  %%  ~Describe the value returned
  \item{gene}{Named chracter vector with mapping from gene ENTREZ id to
    transcript id.}
  \item{exonsNI}{RangedData object containing non overlapping exons.}
  \item{exons}{GRangesList object with original exon coordinates. List
    elements correspond to a known transcripts.}
  \item{txs}{GRanges object with original transcripts annotation.}
  \item{newTxs}{Named list mapping transcript id's to non-overlapping exons.}
  \item{exonmap}{Named list mapping original exon id's to
    non-overlapping exon id's.}


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
\keyword{ ~annotation }
