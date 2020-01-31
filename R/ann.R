ann.reconsensus <- function(ann.data, recode.data, recode.ref.data)
{
	ann.data.crt <- ann.data
	for (i in 1:nrow(ann.data)){
		ann.data.crt$cons_seq[[i]] <- integer()
		for (j in 1:length(ann.data$cons_seq[[i]])){
			cur.prob <- length(findReadsByLoci(recode.data[ann.data$nn_reads_id[[i]]], ann.data$cons_seq[[i]][j])) / length(ann.data$nn_reads_id[[i]])
			if (cur.prob >= 0.1){
				ann.data.crt$cons_seq[[i]] <- c(ann.data.crt$cons_seq[[i]], ann.data$cons_seq[[i]][j])
			}
		}
	}
	ann.data.crt
}

