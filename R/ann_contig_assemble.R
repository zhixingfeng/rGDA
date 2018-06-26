ann_check_nccontig <- function(ann.data)
{
	is_nc <- rep(TRUE, nrow(ann.data))
	for (i in 1:nrow(ann.data)){
		idx <- which(ann.data$start <= ann.data$start[i] & ann.data$end >= ann.data$end[i])
		for (j in idx){
			if (j == i) next
			start_code <- 4*ann.data$start[i]
			end_code <- 4*ann.data$end[i] + 3
			cons_seq_overlap <- ann.data$cons_seq[[j]][ann.data$cons_seq[[j]]>=start_code & ann.data$cons_seq[[j]]<=end_code]
			if (identical(cons_seq_overlap, ann.data$cons_seq[[i]])){
				is_nc[i] <- FALSE
				break
			}
		}
	}
	which(is_nc)
}

ann_contig_to_adjmat <- function(ann.data, min.overlap.prop = 0.5)
{
	adj.mat <- matrix(0, nrow = nrow(ann.data), ncol = nrow(ann.data))
	for (i in 1:nrow(ann.data)){
		overlap.start <- sapply(ann.data$start, function(x,t) max(x,t), t = ann.data$start[i])
		overlap.end <- sapply(ann.data$end, function(x,t) min(x,t), t = ann.data$end[i])
		overlap.len <- overlap.end - overlap.start + 1
		idx <- which(overlap.len >= min.overlap.prop*(ann.data$end[i] - ann.data$start[i] + 1) & ann.data$end >= ann.data$end[i])
		for (j in idx){
			if (j == i) next
			cons_seq_overlap_i <- ann.data$cons_seq[[i]][ann.data$cons_seq[[i]]>=4*overlap.start[j] & ann.data$cons_seq[[i]]<=4*overlap.end[j]+3]
			cons_seq_overlap_j <- ann.data$cons_seq[[j]][ann.data$cons_seq[[j]]>=4*overlap.start[j] & ann.data$cons_seq[[j]]<=4*overlap.end[j]+3]
			if (identical(cons_seq_overlap_i, cons_seq_overlap_j))
				adj.mat[i,j] <- 1
		}
	}
	adj.mat
}


