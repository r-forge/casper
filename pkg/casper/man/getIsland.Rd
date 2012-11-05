\name{getIsland}
\alias{getIsland}
\alias{getIsland,character,missing,annotatedGenome-method}
\alias{getIsland,missing,character,annotatedGenome-method}
\docType{methods}
\title{
getIsland returns the island id associated to a given entrez or
transcript id in an annotatedGenome object.
}
\description{
annotatedGenome objects store information regarding genes and
transcripts. When there's an overlap in exons between several genes,
these genes are grouped into gene islands. getIsland retrieves the
island to which each gene or transcript was assigned.
}
\usage{
getIsland(entrezid, txid, genomeDB)
}
\arguments{
 \item{entrezid}{Character indicating single Entrez identifier. Can be
   left missing and specify txid instead.}
 \item{txid}{Character indicating a single RefSeq transcript identifier. Can be
 left missing and specify entrezid instead}
 \item{genomeDB}{Object of class \code{annotatedGenome}}
}
\section{Methods}{
\describe{
  \item{\code{signature(entrezid='character',txid='missing',genomeDB='annotatedGenome')}}{
    Return island id for given Entrez identifier
  }
  \item{\code{signature(entrezid='missing',txid='character',genomeDB='annotatedGenome')}}{
    Return island id for given transcript identifier (RefSeq)
  }
}}
\value{Character with island identifier}
\examples{
data(hg19DB)
getIsland(entrezid="27",genomeDB=hg19DB)
getIsland(txid="NM_005158",genomeDB=hg19DB)
}
\keyword{manip}
