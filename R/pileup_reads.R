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


