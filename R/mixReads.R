mix.m5 <- function(m5.files, cvg, genome.size, outdir)
{	
	# get average coverage
	m5.cvg <- numeric(length(m5.files))
	for (i in 1:length(m5.files)){
	        print(i)
        	m5.cvg[i] <- get.depth(m5.files[i]) / genome.size
	}
	if (length(cvg) != length(m5.cvg))
		stop("length(cvg) != length(m5.cvg)")
	if (any(cvg>m5.cvg))
		stop("any(cvg>m5.cvg)")
	
	# create output directory
	cmd <- paste('mkdir -p', outdir)		
	print(cmd);system(cmd)

	# sampling 
	for (i in 1:length(m5.files)){
		cur.outfile <- paste(outdir, '/' ,basename(m5.files[i]),'.sample', sep='')
		if (cvg[i]==m5.cvg[i]){
			cmd <- paste('cp', m5.files[i], cur.outfile)
			print(cmd);system(cmd)
		}else{
			m5.reads <- load.m5file(m5.files[i])
			m5.reads.sample <- sample.reads(m5.reads, freq=cvg[i]/m5.cvg[i])
			write.table(m5.reads.sample, file=cur.outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	}
	cmd <- paste('cat ', outdir, '/*.sample > ', outdir, '/pool.m5', sep='')
	print(cmd);system(cmd)
	cmd <- paste('rm ', outdir, '/*.sample', sep='')
	print(cmd);system(cmd)
		
}

sample.reads <- function(m5.reads, freq)
{
	n.reads <- ceiling(freq*nrow(m5.reads))
	idx <- sample(1:nrow(m5.reads), n.reads)
	m5.reads[idx,]
}

