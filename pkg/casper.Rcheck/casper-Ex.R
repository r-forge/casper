pkgname <- "casper"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('casper')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("calcExp")
### * calcExp

flush(stderr()); flush(stdout())

### Name: calcExp
### Title: Calculate expression of gene variants
### Aliases: calcExp
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("genPlot")
### * genPlot

flush(stderr()); flush(stdout())

### Name: genPlot
### Title: Plot exon structure and aligned reads for a given gene
### Aliases: genPlot
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("getDistrs")
### * getDistrs

flush(stderr()); flush(stdout())

### Name: getDistrs
### Title: Compute fragment start and fragment length distributions
### Aliases: getDistrs
### Keywords: ~kwd1 ~kwd2

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




cleanEx()
nameEx("pathCounts")
### * pathCounts

flush(stderr()); flush(stdout())

### Name: pathCounts
### Title: Compute exon path counts
### Aliases: pathCounts

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.





cleanEx()
nameEx("procBam")
### * procBam

flush(stderr()); flush(stdout())

### Name: procBam
### Title: Process SAM/BAM files
### Aliases: procBam
### Keywords: ~paired-end sequencing ~SAM/BAM

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.





cleanEx()
nameEx("procGenome")
### * procGenome

flush(stderr()); flush(stdout())

### Name: procGenome
### Title: Download and format annotations for a given genome
### Aliases: procGenome
### Keywords: ~annotation

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
