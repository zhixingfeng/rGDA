get_nccontig <- function(ann.data, min.prop = 0.5)
{
	is.nc <- rep(TRUE, nrow(ann.data))
	for (i in 1:nrow(ann.data)){
		if (i %% 100 == 0)
			print(i)
		for (j in 1:nrow(ann.data)){
			if (i == j) next
			if (ann.data$start[i] < ann.data$start[j] | ann.data$end[i] > ann.data$end[j])
				next
			overlap.start <- max(ann.data$start[i], ann.data$start[j])
			overlap.end <- min(ann.data$end[i], ann.data$end[j])
			
			cons_seq.i <- ann.data$cons_seq[[i]][ ann.data$cons_seq[[i]]>=4*overlap.start & ann.data$cons_seq[[i]]<=4*overlap.end+3 ]
			cons_seq.j <- ann.data$cons_seq[[j]][ ann.data$cons_seq[[j]]>=4*overlap.start & ann.data$cons_seq[[j]]<=4*overlap.end+3 ]

			n.overlap <- length(intersect(cons_seq.i, cons_seq.j))
			if (n.overlap < min.prop*length(cons_seq.i) | n.overlap < min.prop*length(cons_seq.j)) 
				next
			
			diff.loci.i <- floor(setdiff(cons_seq.i,cons_seq.j)/4)
			diff.loci.j <- floor(setdiff(cons_seq.j,cons_seq.i)/4)
		
			n.diff.loci.i <- length(intersect(diff.loci.i, ann.data$tested_loci[[j]]))
			n.diff.loci.j <- length(intersect(diff.loci.j, ann.data$tested_loci[[i]]))  
			
			if (n.diff.loci.i > 0 | n.diff.loci.j > 0)
				next
			
			is.nc[i] <- FALSE
			break
		}
	}
	idx.nc <- which(is.nc)
	list(ann.data=ann.data[idx.nc,], idx.nc=idx.nc)
}


