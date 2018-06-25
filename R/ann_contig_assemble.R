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
	
}


