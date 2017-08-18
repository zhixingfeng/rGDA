sclust <- function(encode.data, m5.data, subspace, min.count, min.logLR)
{
	if (length(encode.data) != nrow(m5.data))
		stop('size of encode.data and m5.data are not the same.')
	
	rl.count <- rep(0, length(subspace))
	rl.logLR <- rep(0, length(subspace))

	# find reads covering the whole subspace
	locus <- floor(subspace / 4)
	idx <- which(m5.data$tStart<=min(locus) & m5.data$tEnd>=max(locus))
	
	# count reads
	count <- rep(0, 2^length(subspace))
	for (i in 1:length(idx)){
		encode.match <- match(encode.data[[idx[i]]], subspace)
		encode.match <- encode.match[!is.na(encode.match)]
		pattern <- sum(2^(encode.match - 1))
		if (pattern >= 1)
			count[pattern] <- count[pattern] + 1
	}

	# test significance	
	for (i in 1:(length(count)-1)){
		bin <- as.integer(intToBits(i))
		if (sum(bin==1) <= 1 | count[i] < min.count)
			next
		for (j in 1:length(subspace)){
			if (bin[j]==1){
				focal.bit = 2^(j-1)
				cur.logLR <- cal.logLR(count[i], count[focal.bit], count[i-focal.bit], length(idx))
				if (cur.logLR > rl.logLR[j]){
					rl.logLR[j] <- cur.logLR
					rl.count <- count[i]
				}				
			}
		}		
	}
	list(logLR = rl.logLR, count = rl.count)
}

cal.logLR <- function(n11, n10, n01, N)
{
	if (n11/N > (n11+n10)/N * (n11+n01)/N){
		return(cal.logLR.H1(n11, n10, n01, N) - cal.logLR.H0(n11, n10, n01, N))
	}else{
		return(0)
	}
}

cal.logLR.H0 <- function(n11, n10, n01, N)
{
	n00 <- N - n11 - n10 - n01
	n1x <- n11 + n10
	n0x <- n01 + n00
	nx1 <- n11 + n01
	nx0 <- n10 + n00
	
	logLR <- cal.nlogn(n1x) + cal.nlogn(n0x) + cal.nlogn(nx1) + cal.nlogn(nx0) - 2*cal.nlogn(N)				
}

cal.logLR.H1 <- function(n11, n10, n01, N)
{
	n00 <- N - n11 - n10 - n01
	
	logLR <- cal.nlogn(n11) + cal.nlogn(n10) + cal.nlogn(n01) + cal.nlogn(n00) - cal.nlogn(N)
}

cal.nlogn <- function(x)
{
	if (as.integer(x)==0)
		return (0)
	else
		return (x*log(x))
}







