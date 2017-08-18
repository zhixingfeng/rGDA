read.context <- function(context.file)
{
	x <- read.table(context.file, header=FALSE, as.is=TRUE, sep='\t')
	x <- cbind(paste(x[,1],x[,2],sep=','), x)
	names(x) <- c('context', 'context_left', 'context_right', 'A', 'C', 'G', 'T', 'cvg')
	x$context <- as.character(x$context)
	x	
}

