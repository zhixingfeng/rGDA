# detect SNV using single site information and correct context effect
binom.bf <- function(x, n, p0)
{
	a <- 1
	b <- 1
	l1 <- lbeta(x + a, n - x + b) + log(1 - pbeta(p0, x + a, n - x + b) )
	l0 <- x*log(p0) + (n-x)*log(1- p0)
	l1 - l0
}


detect.single <- function(pileup.data, context.effect, min.cvg = 500)
{
	# load context effect 
	context.effect <- context.effect[context.effect$cvg >= min.cvg,]
	context.effect[,4:7] <- context.effect[,4:7] / context.effect$cvg

	# setup result for Bayes factor
	pileup.data.bf <- data.frame(matrix(NaN, nrow(pileup.data), 4))
	
	# find context 
	idx <- match(pileup.data$context, context.effect$context)
        idx.pileup <- which(!is.na(idx))
        idx.context <- idx[!is.na(idx)]
	
	# Bayes binomial test
	if (length(idx.pileup) != length(idx.context))
		stop('length(idx.pileup) != length(idx.context)')

	#rl <- data.frame(matrix(NaN, 4*nrow(pileup.data), 5))
	for (i in 1:length(idx.pileup)){
		if (i %% 1000 == 0)
			print(i)
		cur.pileup.data <- pileup.data[idx.pileup[i], ]
		cur.context.effect <- context.effect[idx.context[i], ]

		# BF for A
		pileup.data.bf[i,1] <- binom.bf(cur.pileup.data$A, cur.pileup.data$cvg, cur.context.effect$A)
		
		# BF for C
                pileup.data.bf[i,2] <- binom.bf(cur.pileup.data$C, cur.pileup.data$cvg, cur.context.effect$C)
	
		# BF for G
                pileup.data.bf[i,3] <- binom.bf(cur.pileup.data$G, cur.pileup.data$cvg, cur.context.effect$G)

		# BF for T
                pileup.data.bf[i,4] <- binom.bf(cur.pileup.data$T, cur.pileup.data$cvg, cur.context.effect$T)
	}
	names(pileup.data.bf) <- c('A_bf', 'C_bf', 'G_bf', 'T_bf')
	
	cbind(pileup.data, pileup.data.bf)
}

#detect.single.print <- function(detect.data, out.file, min.cvg = 20, min.prop = 0.01, min.bf = 50)
#{
#	rl <- data.frame(matrix(NaN, 4*nrow(detect.data), 7))
#	names(rl) <- c('code', 'log_bf', 'prop', 'x', 'n', 'n_link_loci', 'link_loci')
#	rl[,6] <- 0
#	rl[,7] <- -1
#	min.bf.log <- log(min.bf)
#	idx <- which(detect.data$cvg >= min.cvg)
#	for (i in idx){
#		# test A 
#		if (detect.data$A_bf[i] >= min.bf.log & detect.data$A[i]/){
#			t = 4*(i-1) + 1
#			rl$code[t] <- 4*detect.data$locus[i]
#			rl$log_bf <- detect.data$A_bf[i]
#		}
#			
#	}
#}


detect.single.legacy <- function(pileup.data, context.effect, min.cvg = 2000)
{	
	context.effect <- context.effect[context.effect$cvg >= min.cvg,]
	pileup.data.freq <- pileup.data[,6:9]
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
	names(pileup.data)[11:12] <- c('max_lift', 'max_base')
	names(pileup.data)[6:9] <- c('A_lift', 'C_lift', 'G_lift', 'T_lift')
	
	pileup.data <- cbind(pileup.data, pileup.data.freq)
	pileup.data 
		
}





