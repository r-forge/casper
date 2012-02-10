#Remove short inserts
rmShortInserts <- function(bam, isizeMin=100) {
  sel <- abs(bam$isize)>isizeMin
  lapply(bam,function(z) z[sel])
}
