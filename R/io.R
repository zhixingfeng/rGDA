load.encodefile <- function(encodefile)
{
	x.raw <- readLines(encodefile)
	lapply(strsplit(x.raw,'\t'), as.integer )
}

load.m5file <- function(m5file)
{
	x <- read.table(m5file, header=FALSE, as.is=TRUE)
	
	names(x) <- c('qName', 'qLength', 'qStart', 'qEnd', 'qStrand', 'tName', 'tLength', 'tStart', 'tEnd', 'tStrand', 'score', 'numMatch', 'numMismatch', 'numIns', 'numDel', 'mapQV', 'qAlignedSeq', 'matchPattern', 'tAlignedSeq')
	x$tEnd <- x$tEnd - 1
	x$qEnd <- x$qEnd - 1
	x
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

