dist_hamming <- function(encode.1, m5.1, encode.2, m5.2, var.data)
{
	overlap.start <- max(m5.1$tStart, m5.2$tStart)
	overlap.end <- min(m5.1$tEnd, m5.2$tEnd)
	encode.1.overlap <- encode.1[encode.1 >= 4*overlap.start & encode.1 <= 4*overlap.end+3]
	encode.2.overlap <- encode.2[encode.2 >= 4*overlap.start & encode.2 <= 4*overlap.end+3]
	
	n.diff <- length(setdiff(encode.1.overlap, encode.2.overlap)) + length(setdiff(encode.2.overlap, encode.1.overlap))
	n.var <- sum(var.data$locus >= overlap.start & var.data$locus <= overlap.end)

	if (n.var > 0){
		list(dist = n.diff/n.var, n.diff = n.diff, n.var = n.var)
	}else{
		NULL
	}
}




