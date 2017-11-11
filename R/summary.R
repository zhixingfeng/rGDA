summary <- function(sclust.file, out.file, min.logLR, min.count, min.cvg)
{
	if (file.exists(out.file))
		file.remove(out.file)
	sclust.data <- read.table(sclust.file, header=FALSE, as.is=TRUE)
	sclust.data <- sclust.data[sclust.data[,3]>=min.cvg,]
	sclust.data.list <- split(sclust.data[,-1], sclust.data[,1])
	
	for (i in 1:length(sclust.data.list)){
		pattern <- lapply(strsplit(sclust.data.list[[i]][,1],','), as.integer)
		logLR <- lapply(strsplit(sclust.data.list[[i]][,3],','), as.numeric)
		count <- lapply(strsplit(sclust.data.list[[i]][,4],','), as.integer)
		
		pattern.pool <- unlist(pattern)
		logLR.pool <- unlist(logLR)
		count.pool <- unlist(count)
		if (!(length(pattern.pool)==length(logLR.pool) & length(pattern.pool)==length(count.pool)))
			stop("pattern, logLR or count do not match")
		
		pattern.info <- split(as.data.frame(cbind(logLR.pool, count.pool)), pattern.pool)
		logLR.max <- sapply(pattern.info, function(x) max(x[,1]))
		count.max <- sapply(pattern.info, function(x) x[,2][which.max(x[,1])])
		
		idx.sel <- which(logLR.max >= min.logLR & count.max >= min.count)
	
		if (length(idx.sel)==0)
			next
		rl <- paste(as.integer(names(sclust.data.list))[i], '\t', 
			paste(as.integer(names(pattern.info))[idx.sel], collapse=','), ',\t',
			paste(sprintf("%.3f",logLR.max[idx.sel]), collapse=','), ',\t',
			paste(count.max[idx.sel], collapse=','), ',', sep='')
		cat(paste(rl, '\n',sep=''), file=out.file, append=TRUE)	
		
		#rl <- as.data.frame(cbind(as.integer(names(pattern.info))[idx.sel], logLR.max[idx.sel], count.max[idx.sel]))
		
	}	
}



