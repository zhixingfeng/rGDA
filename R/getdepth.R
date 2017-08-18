get.depth <- function(m5file)
{
	cmd <- paste("awk '{sum+=$9-$8} END{print sum}'", m5file)
	as.integer(system(cmd, intern=TRUE))
}

get.depth.samtools <- function(bamfile)
{
	buf <- system(paste('samtools depth',bamfile), intern=TRUE)
	buf.split <- strsplit(buf, '\t')
	data.frame(chromosome=sapply(buf.split,function(x)x[1]), position=as.integer(sapply(buf.split,function(x)x[2])), depth=as.integer(sapply(buf.split,function(x)x[3])))
}



