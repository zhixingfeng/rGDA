encode.data <- load.encodefile('./ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000')
m5.data <- load.m5file('./ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000')
k <- 200
centroid.range <- c(m5.data$tStart[k], m5.data$tEnd[k])
centroid <- encode.data[[k]]
rl <- mat_fac_rank_1_core(encode.data, m5.data, centroid.range, centroid)
rl <- mat_fac_rank_1(encode.data, m5.data, centroid.range, centroid, 10, 1000)

mat_fac_rank_1 <- function(encode.data, m5.data, centroid.range, centroid, min.idx.on = 10, max.iter=1)
{
	old.centroid <- centroid
	n.iter <- 0
	for (i in 1:max.iter){
		rl <- mat_fac_rank_1_core(encode.data, m5.data, centroid.range, old.centroid)
		n.iter <- n.iter + 1
		if (length(rl$idx.on) < min.idx.on)
			break;
			
		if (identical(as.integer(rl$new.centroid), as.integer(old.centroid)))
                        break
                
                old.centroid <- rl$new.centroid
                old.idx.on <- rl$idx.on
		old.idx.off <- rl$idx.off

	}
	rl$n.iter <- n.iter
	rl
}

mat_fac_rank_1_core <- function(encode.data, m5.data, centroid.range, centroid)
{
	###--------- match reads to the centroid --------###
	if (length(encode.data)!=nrow(m5.data)){
		stop('encode.data and m5.data do not match.')
	}
	is.on <- rep(FALSE, length(encode.data))
	for (i in 1:length(encode.data)){
		common.var <- intersect(encode.data[[i]], centroid)
		if (length(common.var) >= ceiling(length(centroid)/2 ))
			is.on[i] <- TRUE
	}
	idx.on <- which(is.on)
	idx.off <- which(!is.on)
	
	###--------- update centroid -------###
	# use the stupid binary matrix to represent reads
	centroid.range.code <- centroid.range
	centroid.range.code[1] <- 4*centroid.range[1]	
	centroid.range.code[2] <- 4*centroid.range[2] + 3
	
	mat.size <- centroid.range.code[2] - centroid.range.code[1] + 1
	
	encode.data.mat.var <- matrix(0, nrow=length(idx.on), ncol=mat.size)
	encode.data.mat.reads <- encode.data.mat.var
	
	for (i in 1:length(idx.on)){
		# fill variants matrix
		idx.mat.var <- intersect(encode.data[[idx.on[i]]],centroid.range.code[1]:centroid.range.code[2]) -
				centroid.range.code[1] + 1
		encode.data.mat.var[i, idx.mat.var] <- 1
		
		# fiil reads cover matrix
		cur.start.code <- 4*m5.data$tStart[idx.on[i]]
		cur.end.code <- 4*m5.data$tEnd[idx.on[i]] + 3
		
		overlap.start.code <- max(cur.start.code - centroid.range.code[1] + 1, 1)
		overlap.end.code <- min(cur.end.code - centroid.range.code[1] + 1, 
					centroid.range.code[2] - centroid.range.code[1] + 1)
		if (overlap.start.code > overlap.end.code)
			stop('overlap.start.code > overlap.end.code')
	
			
		encode.data.mat.reads[i, overlap.start.code:overlap.end.code] <- 1
	}
	# use majority vote to calculate new centroid
	
	encode.data.mat.var.count <- apply(encode.data.mat.var, 2, sum)
	encode.data.mat.reads.count <- apply(encode.data.mat.reads, 2, sum)
	
	if (any(encode.data.mat.var.count > encode.data.mat.reads.count))
		stop('encode.data.mat.var.count > encode.data.mat.reads.count')
	new.centroid <- centroid.range.code[1] + which(encode.data.mat.var.count>=ceiling(encode.data.mat.reads.count/2)
				& encode.data.mat.reads.count>=1) - 1
	list(idx.on=idx.on, idx.off=idx.off, new.centroid=new.centroid)
}



