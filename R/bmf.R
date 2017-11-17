bmf <- function(A, m5.data, init.basis.id, min.overlap=2000)
{
	A.range <- cbind(m5.data$tStart*4, m5.data$tEnd*4 + 3)

	
	y <- A[[init.basis.id]]	
	y.range <- A.range[init.basis.id,]
	x.vec <- integer(length(A))
	
	# update y
	for (i in 1:length(A)){
		overlap.start <- max(min(y.range), min(A.range[i,]))
		overlap.end <- min(max(y.range), max(A.range[i,]))
		if (overlap.end - overlap.start < min.overlap){
			x.vec[i] <- NaN
			next
		}
	
		cur.A <- A[[i]][A[[i]]>=overlap.start & A[[i]]<=overlap.end]
		cur.y <- y[y>=overlap.start & y<=overlap.end]
		H.1 <- length(setdiff(cur.A, cur.y)) + length(setdiff(cur.y, cur.A))
		H.0 <- length(cur.A)
		if (H.1 < H.0){
			x.vec[i] <- 1
		}else{
			x.vec[i] <- 0
		}			
	}

}


