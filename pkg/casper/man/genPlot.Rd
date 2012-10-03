\name{genPlot}
\alias{genPlot}
\title{
  Plot exon structure and aligned reads for a given gene
}
\description{
  Plot exon structure for all transcripts in a given gene and aligned reads
}
\usage{
genPlot(goi, genomeDB, reads, exp)
}
\arguments{
  \item{goi}{
    ENTREZ id of gene of interest
}
  \item{genomeDB}{
    List of annotations produced with the "procGenome" function
}
  \item{reads}{
    RangedData object of aligned reads or fragments 
}
  \item{exp}{
    ExpressionSet object with expression values
}
}
\value{
  \item{exp}{Named integer vector with transcript expression for the
    gene of interest}
}
\author{
  Camille Stephan-Otto Attolini
}
\examples{
##See example in calcExp
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
