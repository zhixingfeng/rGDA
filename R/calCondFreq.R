calCondFreq <- function(m5.data, encode.data, cand.y, cand.x)
{
	cand.yx <- c(cand.y, cand.x)
	idx.yx <- which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = cand.yx))
	idx.x <- which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = cand.x))	
	
	idx.cvg.yx <- which(m5.data[['tStart']] <= floor(min(cand.yx)/4) & m5.data[['tEnd']] >= floor(max(cand.yx)/4))	

	x <- length(idx.yx)
 	n <- length(intersect(idx.x, idx.cvg.yx))
	list(prop = x/n, x = x, n = n)	

}


findReadsByLoci <- function(encode.data, loci)
{
	which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = loci))
}


