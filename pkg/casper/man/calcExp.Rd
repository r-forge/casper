\name{calcExp}
\alias{calcExp}
\title{
  Estimate expression of a known set of gene splicing variants. 
}
\description{
  Estimate expression of gene splicing variants,
  assuming that the set of variants is known.
  When \code{relativeExpr} is set to \code{FALSE}, reads per kilobase
  per million are returned. Otherwise relative expression estimates are returned.
}
\usage{
calcExp(distrs, genomeDB, pc, readLength, islandid, rpkm=TRUE, priorq=2,
priorqGeneExpr=2, citype="none", niter=10^3, burnin=100, mc.cores=1, verbose=FALSE)
}
\arguments{
\item{distrs}{ List of fragment distributions as generated by the \code{getDistrs} function}
\item{genomeDB}{ \code{knownGenome} object containing annotated genome,
  as returned by the \code{procGenome} function. }
\item{pc}{Named vector of exon path counts as returned by \code{pathCounts}}
\item{readLength}{ Read length in bp, e.g. in a paired-end experiment where
  75bp are sequenced on each end one would set \code{readLength=75}.}
\item{islandid}{Name of the gene island to be analyzed. If not specified, all
  gene islands are analyzed.}
\item{rpkm}{Set to \code{FALSE} to return relative expression levels,
  i.e. the proportion of reads generated from each variant per
  gene. These proportions add up to 1 for each gene.
  Set to \code{TRUE} to return reads per kilobase per million (RPKM).}
\item{priorq}{Parameter of the prior distribution on the proportion of reads coming from each variant. The prior is Dirichlet with prior sample size for each variant equal to priorq.
            We recommend \code{priorq=2} for estimation, as it pools the
	    estimated expression away from 0 and 1 and returned lower
	    estimation errors than \code{priorq=1} in our simulated
  experiments.}
\item{priorqGeneExpr}{Parameter for prior distribution on overall gene
  expression. Defaults to 2, which ensures non-zero estimates for all genes}
\item{citype}{Set to \code{"none"} to return no credibility intervals. Set to \code{"asymp"} to return approximate 95\% CIs (obtained via the delta method). Set to \code{"exact"} to obtain exact
            CIs via Monte Carlo simulation. Options \code{"asymp"} and especially \code{"exact"} can increase the computation time substantially.}
\item{niter}{Number of Monte Carlo iterations. Only used when \code{citype=="exact"}.}
\item{burnin}{Number of burnin Monte Carlo iterations. Only used when \code{citype=="exact"}.}
\item{mc.cores}{Number of processors to be used for parallel computation. Can only be used if the package \code{multicore} is available for your system.}
\item{verbose}{Set to \code{TRUE} to display progress information.}
}
\value{
  Expression set with expression estimates.
  \code{featureNames} identify each transcript via
  RefSeq ids, and the \code{featureData} contains further information.
  If \code{citype} was set to a value other than \code{"none"}, the \code{featureData} also contains the 95\% credibility intervals
  (i.e. intervals that contain the true parameter value with 95\% posterior probability).
}
\references{
Rossell D. QASPER: Quantifying Alternative Splicing from Paired End
Reads. Technical report, IRB Barcelona (2010)
}
\author{
  Camille Stephan-Otto Attolini, Manuel Kroiss, David Rossell
}
\examples{
data(K562.r1l1)
data(hg19DB)

#Pre-process
bam0 <- rmShortInserts(K562.r1l1, isizeMin=100)
pbam0 <- procBam(bam0)
head(getReads(pbam0))

#Estimate distributions, get path counts
distrs <- getDistrs(hg19DB,bam=bam0,readLength=75)
pc <- pathCounts(pbam0, DB=hg19DB)

#Get estimates
eset <- calcExp(distrs=distrs, genomeDB=hg19DB, pc=pc, readLength=75, rpkm=FALSE)
head(exprs(eset))
head(fData(eset))

#Add fake sample by permuting and combine
eset2 <- eset[sample(1:nrow(eset),replace=FALSE),]
esetall <- mergeExp(eset,eset2)

#After merge samples are correctly matched
head(exprs(esetall))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% __ONLY ONE__ keyword per line
