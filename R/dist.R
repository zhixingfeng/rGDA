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

sim_jaccard <- function(encode.1, m5.1, encode.2, m5.2, check.ref = FALSE, min.overlap = 500, is.asym = FALSE)
{
	overlap.start <- max(m5.1$tStart, m5.2$tStart)
        overlap.end <- min(m5.1$tEnd, m5.2$tEnd)
        
	if (overlap.end - overlap.start + 1 < min.overlap){
		return(-1)
	}
	
	encode.1.overlap <- encode.1[encode.1 >= 4*overlap.start & encode.1 <= 4*overlap.end+3]
        encode.2.overlap <- encode.2[encode.2 >= 4*overlap.start & encode.2 <= 4*overlap.end+3]

	if (is.asym){
		n.union <- length(union(encode.1, encode.2.overlap))
	}else{
		n.union <- length(union(encode.1.overlap, encode.2.overlap))
	}
	n.intersect <- length(intersect(encode.1.overlap, encode.2.overlap))
	
	
	if (check.ref){
		#if ( n.union == 0 | 2*n.intersect < length(encode.1.overlap) )
		if(n.union == 0 | 2*n.intersect < length(encode.1))
			return(-1)
	}
	#if (n.union == 0 | 2*n.intersect < length(encode.1.overlap) | 2*n.intersect < length(encode.2.overlap) )
	#	return(-1)
	n.intersect / n.union	

}

sim_jaccard_pairwise <- function(cur.encode.data, cur.m5.data, check.ref = FALSE, min.overlap.rate = 0.75)
{
	if (length(cur.encode.data) != nrow(cur.m5.data))
		stop('length(cur.encode.data) != nrow(cur.m5.data)')
	
	dist.mat <- matrix(NaN, length(cur.encode.data), length(cur.encode.data))
	for (i in 1:length(cur.encode.data)){
		if (i %% 100 == 0){
			print (i)
		}
		for (j in 1:length(cur.encode.data)){
			if (i == j) next
			cur.min.overlap <- floor(min.overlap.rate*(cur.m5.data$tEnd[i] - cur.m5.data$tStart[i] + 1))
			dist.mat[i,j] <- sim_jaccard(cur.encode.data[[i]], cur.m5.data[i,], cur.encode.data[[j]], cur.m5.data[j,], check.ref, cur.min.overlap)	
		}
	}
	dist.mat
}


dist_hamming <- function(encode.1, m5.1, encode.2, m5.2, var.data, is.var = FALSE)
{
	overlap.start <- max(m5.1$tStart, m5.2$tStart)
	overlap.end <- min(m5.1$tEnd, m5.2$tEnd)
	encode.1.overlap <- encode.1[encode.1 >= 4*overlap.start & encode.1 <= 4*overlap.end+3]
	encode.2.overlap <- encode.2[encode.2 >= 4*overlap.start & encode.2 <= 4*overlap.end+3]
	
	n.diff <- length(setdiff(encode.1.overlap, encode.2.overlap)) + length(setdiff(encode.2.overlap, encode.1.overlap))
	if (is.var){
		n.var <- sum(var.data$locus >= overlap.start & var.data$locus <= overlap.end)
	}else{
		n.var <- overlap.end - overlap.start + 1
	}

	if (n.var > 0){
		list(dist = n.diff/n.var, n.diff = n.diff, n.var = n.var)
	}else{
		NULL
	}
}




