\name{getNreads}
\alias{getNreads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Function to get number of fragments mapping to each island from an
  object of class 'pathCounts'.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Returns a named numeric vector with number of fragments mapped to each
  island in a pathCounts object. 
}
\usage{
getNreads(pc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pc}{
    %%     ~~Describe \code{pc} here~~
    Object of class 'pathCounts'
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
Named numeric vector. For islands with no path counts (NULL), zero is returned.
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
