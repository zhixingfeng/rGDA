library(dequer)

assemble.core.recode <- function(encode.data, m5.data, centroid, centroid.range, min.overlap=200, min.read.count = 10, is.full.comp = TRUE)
{
	#------------- check input data ------------#
	if (length(encode.data)!=nrow(m5.data)){
                stop('encode.data and m5.data do not match.')
        }

	#------------ match reads to centroids ----------#
        is.on <- rep(FALSE, length(encode.data))
        for (i in 1:length(encode.data)){
		overlap.start <- max(centroid.range[1], m5.data$tStart[i])
		overlap.end <- min(centroid.range[2], m5.data$tEnd[i])
		overlap.len <- overlap.end - overlap.start + 1
		if (overlap.len < min.overlap)
			next

		common.var <- intersect(encode.data[[i]], centroid)
		if (is.full.comp){
			n.centroid.overlap <- length(centroid)
		}else{
                	n.centroid.overlap <- sum(centroid>=4*m5.data$tStart[i] & centroid<=4*m5.data$tEnd[i]+3)
		}
                if (length(common.var)>n.centroid.overlap)
                	stop('length(common.var)>n.centroid.overlap')
                if (length(common.var) >= ceiling(n.centroid.overlap/2) & n.centroid.overlap > 0)
                        is.on[i] <- TRUE
        }
        idx.on <- which(is.on)	

	#------------- update centroid -------------#
	cur.cons <- get_consensus_recode(encode.data[idx.on], m5.data[idx.on,], min.read.count = min.read.count)
	new.centroid <- cur.cons$cons.seq
	
	# restrict new.centroid to centroid.range
	new.centroid <- new.centroid[new.centroid >= 4*centroid.range[1] & new.centroid <= 4*centroid.range[2]+3]

	list(centroid = new.centroid, centroid.range = centroid.range,idx.on = idx.on)
}





