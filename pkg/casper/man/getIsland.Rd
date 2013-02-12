\name{getIsland}
\alias{getIsland}
\alias{getIsland,character,missing,annotatedGenome-method}
\alias{getIsland,missing,character,annotatedGenome-method}
\alias{getChr}
\alias{getChr,character,missing,missing-method}
\alias{getChr,missing,character,missing-method}
\alias{getChr,missing,missing,character-method}
\docType{methods}
\title{
getIsland returns the island id associated to a given entrez or
transcript id in an annotatedGenome object.

getChr indicates the chromosome for a given Entrez, transcript or island id.
}
\description{
annotatedGenome objects store information regarding genes and
transcripts. When there's an overlap in exons between several genes,
these genes are grouped into gene islands. getIsland retrieves the
island to which each gene or transcript was assigned,
while getChr indicates the chromosome.
}
\usage{
getIsland(entrezid, txid, genomeDB)
getChr(entrezid, txid, islandid, genomeDB)
}
\arguments{
 \item{entrezid}{Character indicating single Entrez identifier. Can be
   left missing and specify another identifier instead.}
 \item{txid}{Character indicating a single RefSeq transcript identifier. Can be
   left missing and specify another identifier instead.}
 \item{islandid}{Character indicating the gene island indentifier. Can
   be left missing and specify another identifier instead.}
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
}

\describe{
  \item{\code{signature(entrezid='character',txid='missing',islandid='missing')}}{
    Return chromosome for given Entrez identifier
  }
  \item{\code{signature(entrezid='missing',txid='character',islandid='missing')}}{
    Return chromosome for given transcript identifier (RefSeq)
  }
  \item{\code{signature(entrezid='missing',txid='character',islandid='missing')}}{
    Return chromosome for given island identifier
  }
}

}
\value{Character with island identifier}
\examples{
data(hg19DB)
getIsland(entrezid="27",genomeDB=hg19DB)
getIsland(txid="NM_005158",genomeDB=hg19DB)

getChr(entrezid="27",genomeDB=hg19DB)
getChr(txid="NM_005158",genomeDB=hg19DB)
}
\keyword{manip}
