read.similarity <- function(encode.data, m5.data)
{

}


read.similarity.core <- function(i, j , encode.data, m5.data)
{
	overlap.start <- max(m5.data[i,]$tStart, m5.data[j,]$tStart)
	overlap.end <- min(m5.data[i,]$tEnd, m5.data[j,]$tEnd)
	overlap.start.encode <- overlap.start*4
	overlap.end.encode <- overlap.end*4 + 3
	
	encode.data.1.overlap <- encode.data[[i]][encode.data[[i]]>=overlap.start.encode & encode.data[[i]]<=overlap.end.encode]	
	encode.data.2.overlap <- encode.data[[j]][encode.data[[j]]>=overlap.start.encode & encode.data[[j]]<=overlap.end.encode]	

	encode.data.intersect <- intersect(encode.data.1.overlap, encode.data.2.overlap)
	encode.data.diff.1 <- setdiff(encode.data.1.overlap, encode.data.2.overlap)
	encode.data.diff.2 <- setdiff(encode.data.2.overlap, encode.data.1.overlap)
	
	condprob.diff.1 <- list()
	condprob.diff.2 <- list()
	
	if (length(encode.data.diff.1)>=1){
		for (i in 1:length(encode.data.diff.1))	{
			condprob.diff.1[[i]] <- calCondFreq(m5.data, encode.data, encode.data.diff.1[i], encode.data.intersect)
		}
	}
	if (length(encode.data.diff.2)>=1){
		for (i in 1:length(encode.data.diff.2)) {
                	condprob.diff.2[[i]] <- calCondFreq(m5.data, encode.data, encode.data.diff.2[i], encode.data.intersect)
        	}
	}
	prop.1 <- sapply(condprob.diff.1, function(x) x$prop)
	prop.2 <- sapply(condprob.diff.2, function(x) x$prop)	
	rl <- list()
	rl$condprob.diff.1 <- condprob.diff.1
	rl$condprob.diff.2 <- condprob.diff.2
	rl$encode.data.intersect <- encode.data.intersect
	rl$encode.data.diff.1 <- encode.data.diff.1
	rl$encode.data.diff.2 <- encode.data.diff.2
	rl$prop.1 <- prop.1
	rl$prop.2 <- prop.2
	rl
}


