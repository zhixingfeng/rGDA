library(seqinr)
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

pileup_reads <- function(m5.data)
{
	max.encode <- 4*max(m5.data$tEnd) + 3
	pu.read <- lapply(1:(max.encode+1), function(x) integer(0))
	for (i in 1:nrow(m5.data)){
                if (i %%100 == 0)
                        print(i)
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







