relIsoformExp <- function(eset,geneid='entrezid') {
  if (class(eset)!='ExpressionSet') stop('eset must be of class ExpressionSet')
  geneexpr <- apply(exprs(eset),2,function(z) tapply(z,fData(eset)[,geneid],FUN=sum))
  geneexpr <- geneexpr[as.character(fData(eset)[,geneid]),,drop=FALSE]
  e <- exprs(eset)[,,drop=FALSE]
  sel <- geneexpr!=0
  e[sel,] <- e[sel,,drop=FALSE]/geneexpr[sel,,drop=FALSE]
  exprs(eset) <- e
  eset
}
