\name{simReads}
\alias{simReads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to simulate paired end reads following given read start and
fragment length distributions and gene and variant expressions.
}
\description{
This function generates path counts and bam files with simulated paired
end reads according to given read start distribution, fragment length
distribution and gene and variant expressions.
}
\usage{
simReads(selIslands, geneExpr, rl, n, seed, writeBam, distrs, genomeDB, chrlen=NULL, repSims=FALSE, bamFile=NULL, stranded=FALSE, samtoolsPath=NULL, mc.cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selIslands}{
    Island ID's from the genomeDB object to simulate reads
  }
  \item{geneExpr}{
    Object of class list with the following elements:
    - genes: Numeric vector with relative expressions for each of the simulated islands
      (must add to 1). Names must coincide with 'selIslands'
    - vars: Numeric vector with variant's relative
      expressions. Must add to 1 for each of the simulated genes. Names
      must coincide with transcripts names in the object genomeDB
  }
  \item{rl}{
    Read length
  }
  \item{n}{
    Number of fragments to simulate
  }
  \item{seed}{
    Seed of the random numbers generator
}
  \item{writeBam}{
    Set to 1 to generate bam files with the simulated reads
}
  \item{distrs}{
    Object of class 'readDistrs' with read start and fragment length
    distributions
}
  \item{genomeDB}{
    Object of class 'annotatedGenome' with the genome to genererate
    reads from
}
  \item{repSims}{
    Set to TRUE to return relative read starts and fragment lengths from the simulation
}
  \item{bamFile}{
    Name of the bam file to write reads to. Must end with '.bam' 
}
  \item{stranded}{
    Set to TRUE to preserve gene strand when generating reads. The 'XS'
    tag will be added to reads in the bam file and the returned 'pc'
    object will be stranded
}
  \item{samtoolsPath}{
    Path to the binary file 'samtools'
  }
  \item{mc.cores}{Number of cores to use in function}
  \item{chrlen}{Numeric vector with lengths of chromosomes for reference
    genome. Names must coincide with chromosomes in genomeDB}
}

\value{
  \item{Nsim}{Numerical vector with the number of reads simulated for
    each island.}
  \item{pc}{Object of class 'pathCounts' with simulated path counts}
  \item{sims}{
    Only if 'repSims' is set to TRUE. List with vectors of length 'n'
    with the following elements:
    -'varl': Length of variant for corresponding read
    -'st' Start of fragment relative to variant start (not in genomic
      coordinates)
    -len:Fragment length
    -'strand':Strand of gene for simulated read
  }
}
\author{
Camille Stephan-Otto Attolini
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

data(hg19DB)
data(K562.r1l1)
bam0 <- rmShortInserts(K562.r1l1, isizeMin=100)
distrs <- getDistrs(hg19DB,bam=bam0)

genes <- c('4', '7')
txs <- unlist(lapply(hg19DB@transcripts[genes], names))
vars <- vector(mode='numeric', length=length(txs))
vars[1:3] <- c(1/4,1/4,1/2)
vars[4:7] <- c(0.1, 0.3, 0.5, 0.1)
names(vars) <- txs
geneExpr <- list(genes=c('4'=0.25, '7'=0.75), vars=vars)

simpc <- simReads(selIslands=genes, repSims=TRUE, geneExpr=geneExpr, rl=75,
n=10, seed=1, writeBam=FALSE, distrs=distrs, genomeDB=hg19DB)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
