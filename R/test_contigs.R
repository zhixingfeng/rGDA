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



