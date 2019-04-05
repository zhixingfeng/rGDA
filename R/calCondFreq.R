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


calCondFreqRef <- function(encode.data, encode.ref.data, cand.y, cand.x)
{
        cand.yx <- c(cand.y, cand.x)
        idx.yx <- which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = cand.yx))
        idx.x <- which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = cand.x))
	
	idx.y.A <- which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = 4*floor(cand.y/4) ))
	idx.y.C <- which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = 4*floor(cand.y/4)+1 ))
	idx.y.G <- which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = 4*floor(cand.y/4)+2 ))
	idx.y.T <- which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = 4*floor(cand.y/4)+3 ))
	
	idx.y.A.ref <- which(sapply(encode.ref.data, function(x, t) all(!is.na(match(t,x)) ), t = 4*floor(cand.y/4) ))
        idx.y.C.ref <- which(sapply(encode.ref.data, function(x, t) all(!is.na(match(t,x)) ), t = 4*floor(cand.y/4)+1 ))
        idx.y.G.ref <- which(sapply(encode.ref.data, function(x, t) all(!is.na(match(t,x)) ), t = 4*floor(cand.y/4)+2 ))
        idx.y.T.ref <- which(sapply(encode.ref.data, function(x, t) all(!is.na(match(t,x)) ), t = 4*floor(cand.y/4)+3 ))



	idx.ft <- c(idx.y.A, idx.y.C, idx.y.G, idx.y.T, idx.y.A.ref, idx.y.C.ref, idx.y.G.ref, idx.y.T.ref)


        x <- length(idx.yx)
        n <- length(intersect(idx.x, idx.ft))
        list(prop = x/n, x = x, n = n)

}



