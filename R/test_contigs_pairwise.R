test_contig_pair <- function(ann.data, recode.data, recode.ref.data, i , j)
{
	if (ann.data$start[i] > ann.data$end[j] | ann.data$end[i] < ann.data$start[j])
		return(NULL)

	overlap.start <- max(ann.data$start[i], ann.data$start[j])
	overlap.end <- min(ann.data$end[i], ann.data$end[j])
	overlap.len <- overlap.end - overlap.start + 1

	cons_seq_i <- ann.data$cons_seq[[i]][ ann.data$cons_seq[[i]] >= 4*overlap.start & ann.data$cons_seq[[i]] <= 4*overlap.end+3 ]
	cons_seq_j <- ann.data$cons_seq[[j]][ ann.data$cons_seq[[j]] >= 4*overlap.start & ann.data$cons_seq[[j]] <= 4*overlap.end+3 ]

	fp <- setdiff(cons_seq_i, cons_seq_j)
	fn <- setdiff(cons_seq_j, cons_seq_i)
	diff_var <- c(fp, fn)	

	fn.locus <- floor(fn / 4)

	nr.id <- union(ann.data$nn_reads_id[[i]], ann.data$nn_reads_id[[j]]) + 1
	idx.fp <- findReadsByLoci(recode.data[nr.id], fp)
	idx.fn <- sort(unique( c(findReadsByLoci(recode.ref.data[nr.id], 4*fn.locus), 
				findReadsByLoci(recode.ref.data[nr.id], 4*fn.locus+1),
				findReadsByLoci(recode.ref.data[nr.id], 4*fn.locus+2), 
				findReadsByLoci(recode.ref.data[nr.id], 4*fn.locus+3) ) ))
	cur.count <- length( intersect(idx.fp, idx.fn) )


	diff_var_locus <- floor(diff_var / 4)
	idx.cvg <- nr.id
	for (k in 1:length(diff_var_locus)){
		cur.locus <- diff_var_locus[k]
		cur.idx <- sort(unique( c(findReadsByLoci(recode.data[nr.id], 4*cur.locus),
				findReadsByLoci(recode.data[nr.id], 4*cur.locus+1),
				findReadsByLoci(recode.data[nr.id], 4*cur.locus+2),
				findReadsByLoci(recode.data[nr.id], 4*cur.locus+3),
				findReadsByLoci(recode.ref.data[nr.id], 4*cur.locus),
                                findReadsByLoci(recode.ref.data[nr.id], 4*cur.locus+1),
                                findReadsByLoci(recode.ref.data[nr.id], 4*cur.locus+2),
                                findReadsByLoci(recode.ref.data[nr.id], 4*cur.locus+3) ) ))
		idx.cvg <- intersect(idx.cvg, nr.id[cur.idx])
	}
	cur.cvg <- length(idx.cvg)

	list(count = cur.count, cvg = cur.cvg)
}



