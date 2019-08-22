correct_contig <- function(cur.ann.data)
{
	cons_seq_corrected <- list()
	tested_loci_corrected <- list()
	for (i in 1:length(cur.ann.data$cons_seq)){
		cons_seq_corrected[[i]] <- cur.ann.data$cons_seq[[i]]
		tested_loci_corrected[[i]] <- cur.ann.data$tested_loci[[i]]

		for (j in 1:length(cur.ann.data$cons_seq)){
			if (i == j) next
			locus_j <- floor(cur.ann.data$cons_seq[[j]] / 4)
			cons_seq_j <- cur.ann.data$cons_seq[[j]][ !is.na(match(locus_j, cur.ann.data$tested_loci[[i]])) ]
			if (identical(sort(cur.ann.data$cons_seq[[i]]), sort(cons_seq_j))){
				cons_seq_corrected[[i]] <- sort(union(cons_seq_corrected[[i]], cur.ann.data$cons_seq[[j]]))
				tested_loci_corrected[[i]] <- sort(union(tested_loci_corrected[[i]], floor(cur.ann.data$cons_seq[[j]]/4)))
			}
			
		}
	}
	cur.ann.data$cons_seq <- cons_seq_corrected
	cur.ann.data$tested_loci <- tested_loci_corrected
	cur.ann.data
}







