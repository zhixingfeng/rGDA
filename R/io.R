load.annfile <- function(annfile)
{
	x.raw <- read.table(annfile, header = FALSE, as.is = TRUE, sep = '\t')
	if (ncol(x.raw) == 10 || ncol(x.raw) == 11){
		cons.seq <- lapply(strsplit(as.character(x.raw[,1]),','), as.integer)
                cons.seed <- lapply(strsplit(as.character(x.raw[,8]),','), as.integer)       
                cons.neighbor_id <- lapply(strsplit(as.character(x.raw[,9]),','), as.integer)
                cons.tested_loci <- lapply(strsplit(as.character(x.raw[,10]),','), as.integer)
		
		if (ncol(x.raw) == 11){
			x.raw[,11] <- as.character(x.raw[,11])
			cons.nn_reads_id <- lapply(strsplit(x.raw[,11],','), function(y) if (y[1]==-1) {integer(0)} else {as.integer(y)})
			x <- data.frame(cbind(cons.seq, x.raw[,2], x.raw[,3], x.raw[,4], x.raw[,5], x.raw[,6], x.raw[,7],
                                cons.seed, cons.neighbor_id, cons.tested_loci, cons.nn_reads_id), stringsAsFactors = FALSE)
			names(x) <- c('cons_seq', 'start', 'end', 'contig_count', 'contig_cvg', 'log_bf_null', 'log_bf_ind', 
                                 'seed', 'neighbor_id', 'tested_loci', 'nn_reads_id')
		}else{
			x <- data.frame(cbind(cons.seq, x.raw[,2], x.raw[,3], x.raw[,4], x.raw[,5], x.raw[,6], x.raw[,7],
				cons.seed, cons.neighbor_id, cons.tested_loci), stringsAsFactors = FALSE)		
			names(x) <- c('cons_seq', 'start', 'end', 'contig_count', 'contig_cvg', 'log_bf_null', 'log_bf_ind',
				 'seed', 'neighbor_id', 'tested_loci')
		}

		x$start <- as.integer(x$start)
		x$end <- as.integer(x$end)
		x$contig_count <- as.numeric(x$contig_count)
		x$contig_cvg <- as.numeric(x$contig_cvg)
		x$log_bf_null <- as.numeric(x$log_bf_null)
		x$log_bf_ind <- as.numeric(x$log_bf_ind)
	}
	
	if (ncol(x.raw) == 6){
		cons.seq <- lapply(strsplit(x.raw[,1],','), as.integer)
		cons.seed <- lapply(strsplit(x.raw[,4],','), as.integer)	
		cons.neighbor_id <- lapply(strsplit(x.raw[,5],','), as.integer)
		cons.tested_loci <- lapply(strsplit(x.raw[,6],','), as.integer)
		x <- data.frame(cbind(cons.seq, x.raw[,2], x.raw[,3], cons.seed, cons.neighbor_id, cons.tested_loci), stringsAsFactors = FALSE)
		names(x) <- c('cons_seq', 'start', 'end', 'seed', 'neighbor_id', 'tested_loci')
		x$start <- as.integer(x$start)
		x$end <- as.integer(x$end)
	}
	if (ncol(x.raw) == 3){
		cons.seq <- lapply(strsplit(x.raw[,1],','), as.integer)
		x <- data.frame(cbind(cons.seq, x.raw[,2], x.raw[,3]), stringsAsFactors = FALSE)
		names(x) <- c('cons_seq', 'start', 'end')
		x$start <- as.integer(x$start)
		x$end <- as.integer(x$end)
	}
	return(x)
	
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


load.varfile <- function(var.file)
{
	var.data <- read.table(var.file, as.is = TRUE, sep = '\t')
	names(var.data) <- c('locus', 'base', 'code', 'bf', 'cond_prob', 'joint_prob', 'marginal_prob', 'n_link_loci', 'link_loci')
	var.data
}




