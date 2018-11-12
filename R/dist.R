pair_jaccard <- function(x)
{
	sim_jac <- matrix(0, length(x), length(x))
	for (i in 1:length(x)){
		for (j in 1:length(x)){
			if (j == i) next
			n.common <- length(intersect(x[[i]], x[[j]]))
			n.union <- length(union(x[[i]], x[[j]]))
			sim_jac[i,j] <- n.common / n.union
		}
	}
	sim_jac
}

sim_jaccard <- function(encode.1, m5.1, encode.2, m5.2)
{
	overlap.start <- max(m5.1$tStart, m5.2$tStart)
        overlap.end <- min(m5.1$tEnd, m5.2$tEnd)
        encode.1.overlap <- encode.1[encode.1 >= 4*overlap.start & encode.1 <= 4*overlap.end+3]
        encode.2.overlap <- encode.2[encode.2 >= 4*overlap.start & encode.2 <= 4*overlap.end+3]

	n.union <- length(union(encode.1.overlap, encode.2.overlap))
	n.intersect <- length(intersect(encode.1.overlap, encode.2.overlap))
	
	if (n.union == 0 | 2*n.intersect < length(encode.1.overlap) | 2*n.intersect < length(encode.2.overlap) )
		return(-1)
	n.intersect / n.union	

}

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




