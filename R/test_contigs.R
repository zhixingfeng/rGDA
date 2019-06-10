get_sig_contigs <- function(ann.data, min.log_bf_null = 5, min.log_bf_ind = 5, max.n_loci = 10)
{
	idx.null <- which(ann.data$log_bf_ind == -1000 & ann.data$log_bf_null >= min.log_bf_null)
	idx.ind <- which(ann.data$log_bf_ind >= min.log_bf_ind)
	cons.len <- sapply(ann.data$cons_seq, length)
	idx.max_loci <- which(cons.len >= max.n_loci)
	
	list(idx = union(idx.null, union(idx.ind, idx.max_loci)), idx.null = idx.null, idx.ind = idx.ind, idx.max_loci = idx.max_loci, cons.len = cons.len)	
}


binom_log_bf_bayes <- function(x, n, a_0 = 1.332824, b_0 = 89.04769)
{
	a_1 <- 1
	b_1 <- 1

	log_H1 <- lbeta(x + a_1, n - x + b_1) - lbeta(a_1, b_1)
	log_H0 <- lbeta(x + a_0, n - x + b_0) - lbeta(a_0, b_0)

	log_bf <- log_H1 - log_H0

	log_bf
}

test_contig <- function(cons.seq, pu.var, pu.var.ref)
{
	if (length(cons.seq)<2){
		return(NULL)
	}
	cons.loci <- floor(cons.seq / 4)
	
	pu.cvg <- list()
	for (locus in cons.loci){
		pu.cvg[[locus]] <- integer()
		for (k in 0:3){
			pu.cvg[[locus]] <- union(pu.cvg[[locus]], pu.var[[4*locus+k]])
			pu.cvg[[locus]] <- union(pu.cvg[[locus]], pu.var.ref[[4*locus+k]])
		}
	}

	joint.var <- pu.var[[cons.seq[1]]]
	joint.cvg <- pu.cvg[[cons.loci[1]]]

	for (i in 2:length(cons.seq))
		joint.var <- intersect(joint.var, pu.var[[cons.seq[i]]])
	
	for (i in 2:length(cons.loci))
		joint.cvg <- intersect(joint.cvg, pu.cvg[[cons.loci[i]]])

	prob.obs <- length(joint.var) / length(joint.cvg)

	pu.var.count <- sapply(pu.var[cons.seq], length)
	pu.cvg.count <- sapply(pu.cvg[cons.loci], length)	
	
	prob.exp <- prod(pu.var.count / pu.cvg.count)

	log.bf <- binom.bf(length(joint.var), length(joint.cvg), prob.exp)	
	rr <- prob.obs / prob.exp
		
	list(pu.var.count=pu.var.count, pu.cvg.count=pu.cvg.count, joint.var.count=length(joint.var), joint.cvg.count=length(joint.cvg),
		prob.obs = prob.obs, prob.exp = prob.exp, log.bf = log.bf, rr = rr)
}





