\name{procBam}
\alias{procBam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Process SAM/BAM files
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Process paired-end data stored in SAM or BAM formats and read into
  RangedData objects. The Samtools package is required if files are in
  BAM format. For large files, sequencial processing by chromosome is
  possible in order to minimize RAM needs.
}
\usage{
procBam(bamFileName, samtools, chrom, bam, parallel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bamFileName}{
    %%     ~~Describe \code{bamFileName} here~~
    Absolute path to SAM/BAM file
}
  \item{samtools}{
    %%     ~~Describe \code{samtools} here~~
    Absolute path to the directory where the samtools executable is
}
  \item{chrom}{
    %%     ~~Describe \code{chrom} here~~
    Chromosome to be processed, empty string ("") will process the complete
    genome 
}
  \item{bam}{
    %%     ~~Describe \code{bam} here~~
    0 for SAM format, 1 for BAM format
}
  \item{parallel}{
    %%     ~~Describe \code{parallel} here~~
    Set to TRUE if sequencial processing by chromosome is needed. TRUE
    will overlook the chrom option and process the complete genome.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{reads}{A RangedData object with reads from fragments with both
  ends correctly aligned after splitting them
  by the corresponding CIGAR. Unique identifiers by fragment are
  stored.}
\item{frags}{A RangedData objects with start and end of fragments with
  both reads aligned.}
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
\keyword{ ~paired-end sequencing }
\keyword{ ~SAM/BAM }% __ONLY ONE__ keyword per line
