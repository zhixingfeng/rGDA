findReadsByLoci <- function(encode.data, loci)
{
        which(sapply(encode.data, function(x, t) all(!is.na(match(t,x)) ), t = loci))
}

# calculate distance between the focal locus and link loci in result of dforest
calLinkDist <- function(dforest)
{
	alphabet <- c('A', 'C', 'G', 'T')
	rl <- list()
	rl$focal.locus <- floor(dforest[,1]/4)
	rl$base <- alphabet[(dforest[,1] %% 4)+1]
	rl$condprob <- dforest[,3]
	rl$x <- dforest[,4]
	rl$n <- dforest[,5]
	rl$link.loci <- lapply(strsplit(dforest[,7],','), function(x) floor(as.numeric(x)/4))	
	rl$dist <- mapply(function(x, y) abs(x-y), x=rl$focal.locus, y=rl$link.loci)
	rl$dist.min <- sapply(rl$dist, min)
	rl
}

# get union of intervals
interval_union <- function(x)
{
}

