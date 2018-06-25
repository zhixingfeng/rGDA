dforest_link_loci_olc <- function(dforest.data)
{
	
}

check_link_loci <- function(dforest.data, focal.locus, link.loci)
{
	is.ok <- logical(length(link.loci))
	for (i in 1:length(link.loci)){
		cur.loci <- sort(c(focal.locus, link.loci[-i]))
		idx <- which(dforest.data$focal_locus == link.loci[i])
		is.ok[i] <- any(sapply(dforest.data$link_loci[idx], function(x,t) identical(x,t), t = cur.loci))
	}
}

