correct_reads_jaccard <- function(cur.encode.data, cur.m5.data, check.ref = FALSE, min.overlap.rate = 0.75)
{
	pj <- sim_jaccard_pairwise(cur.encode.data, cur.m5.data, check.ref, min.overlap.rate)
	
	max.common <- list()
        max.range <- list()
	for (i in 1:nrow(pj)){
		max.common[[i]] <- integer()
		max.range[[i]] <- list(start = cur.m5.data$tStart[i], end = cur.m5.data$tEnd[i])
		
		max.jaccard <- max(pj[i,], na.rm = TRUE)
		if (max.jaccard == 0) next
		
		j <- which.max(pj[i,])
                overlap.start <- max(cur.m5.data$tStart[i], cur.m5.data$tStart[j])
                overlap.end <- min(cur.m5.data$tEnd[i], cur.m5.data$tEnd[j])

		if (overlap.start >= overlap.end)
                	stop("overlap.start >= overlap.end")
		
		max.common[[i]] <- intersect(cur.encode.data[[i]], cur.encode.data[[j]])
		max.range[[i]]$start <- overlap.start
                max.range[[i]]$end <- overlap.end	
	}

	list(max.common = max.common, max.range = max.range, pj = pj)
}


correct_reads_maxlen <- function(cur.encode.data, cur.m5.data)
{
	max.common <- list()
	max.range <- list()
	for (i in 1:length(cur.encode.data)){
		if (i %% 100 == 0) print (i)
		max.common[[i]] <- integer()
		max.range[[i]] <- list(start = cur.m5.data$tStart[i], end = cur.m5.data$tEnd[i])
		for (j in 1:length(cur.encode.data)){
			if (i == j) next
			cur.common <- intersect(cur.encode.data[[i]], cur.encode.data[[j]])
			overlap.start <- max(cur.m5.data$tStart[i], cur.m5.data$tStart[j])
			overlap.end <- min(cur.m5.data$tEnd[i], cur.m5.data$tEnd[j])
			if (length(cur.common) > length(max.common[[i]])){
				max.common[[i]] <- cur.common
				if (overlap.start >= overlap.end)	
					stop("overlap.start >= overlap.end")	
				max.range[[i]]$start <- overlap.start
				max.range[[i]]$end <- overlap.end
			}
		}
	}
	list(max.common = max.common, max.range = max.range)
}

