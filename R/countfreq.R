countfreq <- function(encode.data, m5.data, subspace, pattern)
{
	subspace.range <- floor(range(subspace)/4)	
		
	idx.cover <- which(m5.data$tStart <= subspace.range[1] & m5.data$tEnd >= subspace.range[2])
	
	pattern.bin <- rev(as.binary(pattern))
	idx.var <- which(sapply(encode.data[idx.cover], function(x,t,p) identical(intersect(x,p),t), t=subspace[which(pattern.bin)], p=subspace))	
	
	c(length(idx.var), length(idx.cover))
}


