\name{procGenome}
\alias{procGenome}
\alias{createDenovoGenome}
\title{
  Create an annotatedGenome object
}
\description{
  \code{procGenome} downloads UCSC annotations for a given genome.
  \code{createDenovoGenome} creates a de novo annotated genome by
  combining UCSC annotations and observed RNA-seq data.
}
\usage{
procGenome(genome, mc.cores = 1)

createDenovoGenome(reads, DB, readLen, stranded=FALSE, minLinks=2,
maxLinkDist=1e+05, maxDist=1000, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genome}{Character indicating genome version (e.g. "hg19", "dm3")}
  \item{mc.cores}{Number of cores to use in parallel processing
    (multicore package required)}
  \item{DB}{\code{annotatedGenome} object, as returned by
    \code{procGenome}}
  \item{readLen}{Read length, e.g. 75 for paired-end experiments where
    75bp from each end were sequenced}
  \item{stranded}{Set to \code{TRUE} if a stranded RNA-seq experiment
    was performed, e.g. if we know the strand of the transcript that
    generated each read. Set to \code{FALSE} for non-stranded
    experiments.}
  \item{minLinks}{Minimum number of reads joining two
    exons to merge their corresponding genes.}
  \item{maxLinkDist}{Maximum distance between two exons to merge
    their correspondin genes. A value of \code{0} disables this option.}
  \item{maxDist}{Maximum distance between two exons with reads joining
    them to merge their corresponding genes.}
  \item{reads}{Processed reads stored in a \code{RangedData}, as
    returned by \code{procBam}}
}
\details{
  These functions create the annotation objects that are needed for subsequent
  functions. Typically these objects are created only once for a set of
  samples. For \code{createDenovoGenome} it is convenient to create a
  single genome by combining information from all samples, as this
  increases the power to detect new exons and
  guarantees that the same exons and gene islands are used across samples.
}
\value{
  Object of class \code{annotatedGenome}.
}
\author{
  %%  ~~who you are~~
  Camille Stephan-Otto Attolini
}
\seealso{
  annotatedGenome-class for a description of the class.
  To obtain the island id for an entrez gene id or a transcript id see help(getIsland).
  To extract the transcripts for a entrez or island id see help(transcripts).
}
\examples{
## hg19DB <- procGenome("hg19")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~annotation }
