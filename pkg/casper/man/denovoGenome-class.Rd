\name{denovoGenome-class}
\Rdversion{1.1}
\docType{class}
\alias{denovoGenome-class}
\alias{show,denovoGenome-method}

\title{Class "denovoGenome"}
\description{
  Stores de novo annotation for a genome,
  which was typically obtained by combining a known genome with observed
  RNA-seq data.
}
\section{Objects from the Class}{
Objects are typically created with a call to \code{createDenovoGenome}, which
downloads genome annotation from UCSC.
}
\section{Slots}{
  \describe{
    \item{genes}{Each element in the list corresponds to a gene island
      and contains an \code{IRanges} object. It indicates the
      start/end/exon id of each exon contained in the island}
    \item{transcripts}{Each element in the list corresponds to a gene island. It indicates the exons contained in each known variant.}
    \item{exon2gene}{\code{data.frame} indicating the chromosome, start and end of each exon, and its corresponding gene island.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "denovoGenome")}: Displays general information about the object. }
	 }
}
\author{
Camille Stephan-Otto Attolini
}
\seealso{
  \code{\link{knownGenome-class}} stores known annotated genomes, e.g. as contained in the UCSC repository.
}
\examples{
showClass("denovoGenome")
}
\keyword{classes}
