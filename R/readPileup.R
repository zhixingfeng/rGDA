read.pileup <- function(pileupfile)
{
	x <- read.table(pileupfile, as.is=TRUE, header=FALSE, sep='\t')	
	x <- cbind(x[,1], paste(x[,2],x[,3],sep=','), x[,2:9])
	names(x) <- c('locus', 'context', 'context_left', 'context_right', 'ref', 'A', 'C', 'G', 'T', 'cvg')
	x$context <- as.character(x$context)
	x
}



