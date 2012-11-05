\name{procBam-class}
\Rdversion{1.1}
\docType{class}
\alias{procBam-class}
\alias{show,procBam-method}

\title{Class "procBam"}
\description{
  Stores processed bam files in a \code{RangedData} format. Each read is
  split into disjoint ranges according to its cigar code.
}
\section{Objects from the Class}{
  Objects are created with a call to \code{procBam}.
}
\section{Slots}{
  \describe{
    \item{pbam}{\code{RangedData} indicating chromosome, start and end
      of each disjoint range. The pair id and read id within the pair
      are also stored.

    }
    \item{stranded}{Logical variable. \code{TRUE} indicates that the
      reads were obtained from and RNA-seq experiment where strand
      information was preserved.
    }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "procBam")}: Displays general information about the object. }
	 }
}
\author{
Camille Stephan-Otto Attolini, David Rossell
}
\examples{
showClass("procBam")
}
\keyword{classes}
