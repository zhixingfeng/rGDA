library(seqinr)

pileup_var_count_recode <- function(recode.data, recode.ref.data, var.data)
{
	pu_var <- rep(0, 4*max(var.data$locus)+3)
	pu_var_ref <- rep(0, 4*max(var.data$locus)+3)
	
	pu_var[1:max(unlist(recode.data))] <- pileup_var_count(recode.data)
	pu_var_ref[1:max(unlist(recode.ref.data))] <- pileup_var_count(recode.ref.data)
	
	cvg <- pu_var[4*var.data$locus] + pu_var[4*var.data$locus + 1] + pu_var[4*var.data$locus + 2] + pu_var[4*var.data$locus + 3]
	cvg <- cvg + pu_var_ref[4*var.data$locus] + pu_var_ref[4*var.data$locus + 1] + pu_var_ref[4*var.data$locus + 2] + pu_var_ref[4*var.data$locus + 3]
	
	prop.A <- pu_var[4*var.data$locus] / cvg
	prop.C <- pu_var[4*var.data$locus + 1] / cvg 	
	prop.G <- pu_var[4*var.data$locus + 2] / cvg
	prop.T <- pu_var[4*var.data$locus + 3] / cvg
	
	rl <- data.frame(cbind(var.data$locus, prop.A, prop.C, prop.G, prop.T, 
			pu_var[4*var.data$locus], pu_var[4*var.data$locus+1], pu_var[4*var.data$locus+2], 
			pu_var[4*var.data$locus+3], cvg), stringsAsFactors=FALSE)
	
}

pileup_var_count <- function(encode.data)
{	
	max.encode <- max(unlist(encode.data))
	var.count <- rep(0,max.encode)
	for (i in 1:length(encode.data)){
		var.count[encode.data[[i]]]  <- var.count[encode.data[[i]]] + 1
	}
	var.count
}


pileup_reads_count <- function(m5.data){
	max.encode <- 4*max(m5.data$tEnd) + 3
	read.count <- rep(0,max.encode + 1)
	for (i in 1:nrow(m5.data)){
		if (i %%1000 == 0)
			print(i)
		# remove insertions
		idx <- which(s2c(m5.data$tAlignedSeq[i]) != '-')
		t.align <- s2c(m5.data$tAlignedSeq[i])[idx]
		q.align <- s2c(m5.data$qAlignedSeq[i])[idx]

		loci.covered <- m5.data$tStart[i] + which(q.align != '-') - 1
		idx.encode <- c(4*loci.covered, 4*loci.covered+1, 4*loci.covered+2, 4*loci.covered+3)
		read.count[idx.encode + 1] <- read.count[idx.encode + 1] + 1
	}
	print(i)
	read.count
}

pileup_reads <- function(m5.data, min.len = 500)
{
	max.encode <- 4*max(m5.data$tEnd) + 3
	pu.read <- lapply(1:(max.encode+1), function(x) integer(0))
	for (i in 1:nrow(m5.data)){
                if (i %%100 == 0)
                        print(i)
                # discard short reads
		if (m5.data$tEnd[i] - m5.data$tStart[i] +1 < min.len)
			next
		# remove insertions
                idx <- which(s2c(m5.data$tAlignedSeq[i]) != '-')
                t.align <- s2c(m5.data$tAlignedSeq[i])[idx]
                q.align <- s2c(m5.data$qAlignedSeq[i])[idx]

                loci.covered <- m5.data$tStart[i] + which(q.align != '-') - 1
                for (j in 1:length(loci.covered)){
			pu.read[[4*loci.covered[j] + 1]] <- c(pu.read[[4*loci.covered[j] + 1]], i-1)
			pu.read[[4*loci.covered[j] + 1 + 1]] <- c(pu.read[[4*loci.covered[j] + 1 + 1]], i-1)
			pu.read[[4*loci.covered[j] + 2 + 1]] <- c(pu.read[[4*loci.covered[j] + 2 + 1]], i-1)
			pu.read[[4*loci.covered[j] + 3 + 1]] <- c(pu.read[[4*loci.covered[j] + 3 + 1]], i-1)
		}
	}
        print(i)
	pu.read
}

pileup_var <- function(encode.data, max.encode)
{
	var.count <- lapply(1:(max.encode+1), function(x) integer(0))
	for (i in 1:length(encode.data)){
		if (length(encode.data[[i]]) == 0)
			next
		for (j in 1:length(encode.data[[i]])){
			var.count[[encode.data[[i]][j] + 1]] <- c(var.count[[encode.data[[i]][j] + 1]], i - 1)
		}
	}
	var.count
}


print_pileup <- function(pu, out.file )
{
	pu.string <- sapply(pu, paste, collapse = ',')
	writeLines(pu.string, out.file)
}

filter_pileup_var <- function(pu.var, pu.reads)
{
	if (length(pu.var) != length(pu.reads))
		stop('length(pu.var) != length(pu.reads)')
	pu.var.ft <- list()
	for (i in 1:length(pu.var)){
		pu.var.ft[[i]] <- intersect(pu.var[[i]], pu.reads[[i]])
	}
	pu.var.ft
}




