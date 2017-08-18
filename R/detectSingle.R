# detect SNV using single site information and correct context effect
detect.single <- function(pileup.data, context.effect, min.cvg = 2000)
{
	context.effect <- context.effect[context.effect$cvg >= min.cvg,]
	pileup.data[,6:9] <- pileup.data[,6:9] / pileup.data$cvg
	context.effect[,4:7] <- context.effect[,4:7] / context.effect$cvg

	idx <- match(pileup.data$context, context.effect$context)	
	idx.pileup <- which(!is.na(idx))
	idx.context <- idx[!is.na(idx)]

	
	pileup.data[idx.pileup, 6:9] <- pileup.data[idx.pileup, 6:9] / context.effect[idx.context, 4:7]
	pileup.data[is.na(idx),6:9] <- NaN
	
	max.lift <- apply(pileup.data[, 6:9], 1, function(x) {if(all(is.na(x))) return(NaN); max(x,na.rm=TRUE)} )
	max.base <- c('A','C','G','T')[apply(pileup.data[, 6:9],
			1, function(x) {if(all(is.na(x))) return(NaN); which.max(x)} )]
	pileup.data <- cbind(pileup.data, max.lift, max.base)
	names(pileup.data)[11:12] <- c('max.lift', 'max.base')
	
	
	pileup.data
		
}

