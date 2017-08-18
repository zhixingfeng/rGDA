# detect SNV using multiple sites information and correct context effect 
detect.multiple <- function(dforest.max.data, pileup.data, context.effect, min.cvg = 2000, alpha=0.005)
{
		
	# filter context effect by minimal coverage
	context.effect <- context.effect[context.effect$cvg >= min.cvg,]
	context.effect[,4:7] <- context.effect[,4:7] / context.effect$cvg

	# parse dforest.max.data
	link.dist <- calLinkDist(dforest.max.data)
	
	# find context in pile.data(genome coordinate in dforest.max.data is 0-based and in pileup.data is 1-based!!)
	context <- pileup.data$context[match(link.dist$focal.locus+1, pileup.data$locus)]	

	idx.context <- match(context, context.effect$context)	
					
	idx.A <- which(link.dist$base=='A')
	idx.C <- which(link.dist$base=='C')
	idx.G <- which(link.dist$base=='G')
	idx.T <- which(link.dist$base=='T')
	
	lift <- rep(NA, length(link.dist))	
	lift[idx.A] <- link.dist$condprob[idx.A] / (context.effect$A[idx.context[idx.A]] + alpha)
	lift[idx.C] <- link.dist$condprob[idx.C] / (context.effect$C[idx.context[idx.C]] + alpha)
	lift[idx.G] <- link.dist$condprob[idx.G] / (context.effect$G[idx.context[idx.G]] + alpha)
	lift[idx.T] <- link.dist$condprob[idx.T] / (context.effect$T[idx.context[idx.T]] + alpha)
	
	condprob.null <- rep(NA, length(link.dist))
	condprob.null[idx.A] <- context.effect$A[idx.context[idx.A]]
	condprob.null[idx.C] <- context.effect$C[idx.context[idx.C]]
	condprob.null[idx.G] <- context.effect$G[idx.context[idx.G]]
	condprob.null[idx.T] <- context.effect$T[idx.context[idx.T]]	

	link.dist$lift <- lift
	link.dist$context <- context
	link.dist$condprob.null <- condprob.null
	link.dist
}

