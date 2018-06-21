load.annfile <- function(annfile)
{
	x.raw <- read.table(annfile, header = FALSE, as.is = TRUE, sep = '\t')
	cons.seq <- lapply(strsplit(x.raw[,1],','), as.integer)
	cons.seed <- lapply(strsplit(x.raw[,2],','), as.integer)	
	cons.loci <- lapply(strsplit(x.raw[,3],','), as.integer)
	cons.prop <- lapply(strsplit(x.raw[,4],','), as.numeric)
	cons.pu_var_count <- lapply(strsplit(x.raw[,5],','), as.integer)
	cons.pu_read_count <- lapply(strsplit(x.raw[,6],','), as.integer)
	x <- data.frame(cbind(cons.seq, cons.seed, cons.loci, cons.prop, cons.pu_var_count, cons.pu_read_count), stringsAsFactors = FALSE)
	names(x) <- c('cons_seq', 'seed', 'loci', 'prop', 'pu_var_count', 'pu_read_count')
	x
}

load.jaccardfile <- function(jaccardfile)
{
	x <- read.table(jaccardfile, header=FALSE, as.is=TRUE, sep=',')	
	names(x) <- c('id_1', 'id_2', 'jaccard', 'n_common', 'n_union', 'overlap', 'start_1', 'end_1', 'start_2', 'end_2')
	x
}

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
	m5.data$tEnd <- m5.data$tEnd + 1
	m5.data$qEnd <- m5.data$qEnd + 1
	write.table(m5.data, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ' ', file = out.file )
}

load.m4file <- function(m4file){
	x <- read.table(m4file, header=FALSE, as.is=TRUE)
	names(x) <- c('qName', 'tName', 'score', 'percentSimilarity', 'qStrand', 'qStart', 'qEnd', 'qLength', 'tStrand', 'tStart', 'tEnd', 'tLength', 'mapQV')
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

load.testfile <- function(test.file)
{
	test.data.raw <- read.table(test.file, header = FALSE, sep = '\t', as.is = TRUE)
	test.data <- test.data.raw[,1:3]
	names(test.data) <- c('read_id', 'start', 'end')

	test.data$common <- lapply(strsplit(test.data.raw[,4],','), as.integer) 
	test.data$diff <- lapply(strsplit(test.data.raw[,5],','), as.integer)

	test.data$common_prop <- lapply(strsplit(test.data.raw[,6],','), as.numeric)
	test.data$diff_prop <- lapply(strsplit(test.data.raw[,7],','), as.numeric)
	test.data
}


