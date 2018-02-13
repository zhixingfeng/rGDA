check_contained_reads <- function(encode.data, m5.data, min.overlap = 0, rm_empty_centroid = TRUE)
{
	###--------- check if encode.data and m5.data match --------###
        if (length(encode.data)!=nrow(m5.data)){
                stop('encode.data and m5.data do not match.')
        }

	###--------- remove reads contained in another read ------###
	is.contained <- rep(FALSE, length(encode.data))
	for (i in 1:length(encode.data)){
                if (i %% 100 == 0) cat(i, '\r')
		if (length(encode.data[[i]]) == 0){
			if (rm_empty_centroid)
				is.contained[i] <- TRUE
			next
		}
		for (j in 1:length(encode.data)){
                        if (i==j) next
                        if (m5.data$tStart[i] < m5.data$tStart[j] | m5.data$tEnd[i] > m5.data$tEnd[j])
                                next
			if (m5.data$tStart[i] == m5.data$tStart[j] & m5.data$tEnd[i] == m5.data$tEnd[j])
				next
			if (m5.data$tEnd[i] - m5.data$tStart[i] + 1 <  min.overlap)
				next
                        
			n.common.var <- length(intersect(encode.data[[i]], encode.data[[j]]))
                        n.var.overlap <- sum(encode.data[[j]]>=4*m5.data$tStart[i] & encode.data[[j]]<=4*m5.data$tEnd[i]+3)
			if (n.var.overlap == 0)
				next
				
                        if (n.common.var >= ceiling(n.var.overlap/2)){
                        	is.contained[i] <- TRUE
				break
			}
                }

	}
	cat(i, '\n')
	which(!is.contained)
}

