\name{knownGenome-class}
\Rdversion{1.1}
\docType{class}
\alias{knownGenome-class}
\alias{show,knownGenome-method}

\title{Class "knownGenome"}
\description{
  Stores genome annotation for a known genome.
  This class is useful to store known genes, transcripts and exons that
  was downloaded from some public database, e.g. UCSC.
}
\section{Objects from the Class}{
Objects are typically created with a call to \code{procGenome}, which
downloads genome annotation from UCSC.
}
\section{Slots}{
  \describe{
    \item{gene}{List of gene ids}
    \item{exonsNI}{\code{RangedData} indicating start and end of each
      exon}
    \item{exons}{\code{GRangesList} indicating the exons for each known
      transcript}
    \item{txs}{\code{GRanges} object indicating for each transcript the 
      chromosome, strand, transcript id and name. It also stores the
      corresponding gene and exon ids.}
    \item{newTxs}{Each element in the list corresponds to a transcript,
      and contains its exon ids}
    \item{exonmap}{Each element in the list corresponds to an exon}
    \item{gene2ex}{Each element in the list corresponds to a gene island
    and indicates its exon ids}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "knownGenome")}: Displays general information about the object. }
  }
}
\author{
Camille Stephan-Otto Attolini
}
\seealso{
  \code{\link{denovoGenome-class}} stores de novo annotated genomes,
  which combine known annotation with new exons/transcripts inferred
  from the observed data.
}
\examples{
showClass("knownGenome")
}
\keyword{classes}
