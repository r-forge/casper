\name{wrapKnown}
\alias{wrapKnown}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Run all necessary steps to analyze a bam file with the casper pipeline.
}
\description{
  Function to analyze a bam file to generate expression estimates,
  read start and fragment length distributions, path counts and
  processed reads.
}
\usage{
wrapKnown(bamFile, verbose = FALSE, seed = 1, mc.cores.int = 1, mc.cores
= 1, genomeDB, readLength, rpkm = TRUE, priorq = 2, priorqGeneExpr = 2,
citype = "none", niter = 10^3, burnin = 100, keep.pbam = FALSE)
}
\arguments{
  \item{bamFile}{
    Bam file with the sample to analyze. Index of bam file must be in
    the same directory.
}
  \item{verbose}{
    Set to \code{TRUE} to display progress information.
}
  \item{seed}{
    Set seed of random number generator.
}
  \item{mc.cores.int}{
    Number of cores to use when loading bam files. This is a memory intensive
    step, therefore number of cores must be chosen according to available RAM
    memory.
}
  \item{mc.cores}{
    Number of cores to use in expression estimation.
}
  \item{genomeDB}{
    \code{knownGenome} object containing annotated genome,
    as returned by the \code{procGenome} function.
}
  \item{readLength}{
    Read length in bp, e.g. in a paired-end experiment where
    75bp are sequenced on each end one would set \code{readLength=75}.
}
  \item{rpkm}{
    Set to \code{FALSE} to return relative expression levels,
    i.e. the proportion of reads generated from each variant per
    gene. These proportions add up to 1 for each gene.
    Set to \code{TRUE} to return reads per kilobase per million (RPKM).
}
  \item{priorq}{
    Parameter of the prior distribution on the proportion of reads
    coming from each variant. The prior is Dirichlet with prior sample
    size for each variant equal to priorq.
    We recommend \code{priorq=2} for estimation, as it pools the
    estimated expression away from 0 and 1 and returned lower
    estimation errors than \code{priorq=1} in our simulated
    experiments.
}
  \item{priorqGeneExpr}{
    Parameter for prior distribution on overall gene
    expression. Defaults to 2, which ensures non-zero estimates for
    all genes
}
  \item{citype}{
    Set to \code{"none"} to return no credibility intervals. Set to
    \code{"asymp"} to return approximate 95\% CIs (obtained via the
    delta method). Set to \code{"exact"} to obtain exact
    CIs via Monte Carlo simulation. Options \code{"asymp"}
    and especially \code{"exact"} can increase the
    computation time substantially.
}
  \item{niter}{
    Number of Monte Carlo iterations. Only used when \code{citype=="exact"}.
}
  \item{burnin}{
    Number of burnin Monte Carlo iterations. Only used when \code{citype=="exact"}.
  }
  \item{keep.pbam}{Set to \code{TRUE} to save processed bam object, as
    returned by \code{procBam}. This object can require substantial memory
    during execution and disk storage upon saving and is not needed
    for a default analysis.}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  Function executes the functions procBam, getDistrs and pathCounts
  for each chromosome. Once finished, distributions are merged and
  expressions computed. Distributions, processed reads, path counts and
  expression estimates are returned.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{distr}{Object of class \code{readDistrs}}
\item{pbam}{List of objects of class \code{procBam} with one element
  per chromosome}
\item{pc}{Object of class \code{pathCounts}}
\item{exp}{Object of class \code{expressionSet}}
}
\author{
  %%  ~~who you are~~
  Camille Stephan-Otto Attolini
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## genDB<-makeTranscriptDbFromUCSC(genome="hg19", tablename="refGene")
## hg19DB <- procGenome(genDB, "hg19")
  bamFile="/path_to_bam/sorted.bam"
## ans <- wrapKnown(bamFile=bamFile, mc.cores.int=4, mc.cores=3, genomeDB=hg19DB, readLength=101)
##  names(ans)
##  head(exprs(ans\$exp))
}
