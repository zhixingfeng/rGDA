load.encodefile <- function(encodefile)
{
	x.raw <- readLines(encodefile)
	lapply(strsplit(x.raw,'\t'), as.integer )
}

save.encodefile <- function(encode.data, out.file)
{
	encode.data.txt <- sapply(encode.data, function(x) paste(x,collapse='\t'))
	writeLines(encode.data.txt, out.file)
}

load.m5file <- function(m5file)
{
	x <- read.table(m5file, header=FALSE, as.is=TRUE)
	
	names(x) <- c('qName', 'qLength', 'qStart', 'qEnd', 'qStrand', 'tName', 'tLength', 'tStart', 'tEnd', 'tStrand', 'score', 'numMatch', 'numMismatch', 'numIns', 'numDel', 'mapQV', 'qAlignedSeq', 'matchPattern', 'tAlignedSeq')
	x$tEnd <- x$tEnd - 1
	x$qEnd <- x$qEnd - 1
	x
}

save.m5file <- function(m5.data, out.file)
{
	write.table(m5.data, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ' ', file = out.file )
}

load.cmpreadsfile <- function(cmpreadsfile, is.full=FALSE)
{
	x.raw <- readLines(cmpreadsfile)
	lapply(strsplit(x.raw,','), as.integer )
}

load.output <- function(outputfile)
{
	x.raw <- readLines(outputfile)
	x <- strsplit(x.raw, '\t')
	x
}

load.dist <- function(dist.file, min.overlap = 100){
	dist.raw <- read.table(dist.file, sep=',', header=FALSE, as.is=TRUE)
	n.reads <- max(dist.raw[,1] + 1, dist.raw[,2] + 1)	
	dist.mat <- matrix(1, nrow=n.reads, ncol=n.reads)
	for (i in 1:nrow(dist.raw)){
		if (i %% 10000 == 0)
        	        cat(i, '\n')
	        idx.row <- dist.raw[i,1] + 1
        	idx.col <- dist.raw[i,2] + 1
		if (dist.raw[i,5] >= min.overlap)
			dist.mat[idx.row, idx.col] <- dist.raw[i,3]
	}
	diag(dist.mat) <- 0
	dist.mat
}


