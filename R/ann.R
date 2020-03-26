check.ann <- function(ann.data)
{
	
	tested.loci.min <- sapply(ann.data$tested_loci, min)
	tested.loci.max <- sapply(ann.data$tested_loci, max)

	is.start.ok = all(tested.loci.min >= ann.data$start)
	is.end.ok <- all(tested.loci.max <= ann.data$end)	
	diff.start = tested.loci.min - ann.data$start
	diff.end <- ann.data$end - tested.loci.max
	list(is.start.ok = is.start.ok, is.end.ok = is.end.ok, diff.start = diff.start, diff.end = diff.end)
} 

ann.reconsensus <- function(ann.data, recode.data, recode.ref.data,  var.data, min.prop = 0.5, min.cvg = 10)
{
	ann.data.crt <- ann.data
        for (i in 1:nrow(ann.data)){
		cur.idx <- ann.data$nn_reads_id[[i]] + 1
		cur.pu <- pileup_var_count_recode(recode.data[cur.idx], recode.ref.data[cur.idx],  var.data)
		cur.pu.trim <- cur.pu[cur.pu$locus >= ann.data$start[i] & cur.pu$locus <= ann.data$end[i], ]
		ann.data.crt$cons_seq[[i]] <- integer()
		for (j in 1:nrow(cur.pu.trim)){
			cur.cvg <- cur.pu.trim$cvg[j]
			if (cur.cvg < min.cvg) next
			
			cur.max.prop <- max(cur.pu.trim[j, 2:5])
			cur.base <- which.max(cur.pu.trim[j, 2:5]) - 1
			if (cur.max.prop > min.prop){
				cur.code <- 4*cur.pu.trim$locus[j] + cur.base
				ann.data.crt$cons_seq[[i]] <- c(ann.data.crt$cons_seq[[i]], cur.code)
			}
		}
		names(ann.data.crt$cons_seq[[i]]) <- NULL
	}
	ann.data.crt
}


ann.reconsensus.bak <- function(ann.data, recode.data, recode.ref.data)
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

