#find_contained_reads <- function(encode.data, m5.data, min.overlap.var = 1, min.overlap.len=200)
#{
	###--------- check if encode.data and m5.data match --------###
#        if (length(encode.data)!=nrow(m5.data)){
#                stop('encode.data and m5.data do not match.')
#        }

        ###--------- remove reads contained in another read ------###
#        read.pair.mat <- matrix(FALSE, length(encode.data), length(encode.data))
#        for (i in 1:length(encode.data)){
#                print(i)
#                for (j in 1:length(encode.data)){
                        # don't compare reads to themselves
#			if (i==j) next
			# don't compare non-overlaping reads
#                       if (m5.data$tStart[i] > m5.data$tEnd[j] | m5.data$tEnd[i] < m5.data$tStart[j])
#                                next
#			# get overlap region
#			overlap.start <- max(m5.data$tStart[i], m5.data$tStart[j])	
#			overlap.end <- min(m5.data$tEnd[i], m5.data$tEnd[j])
#			overlap.len <- overlap.end - overlap.start + 1
#			if (overlap.len < min.overlap.len)
#				next
			
			# get number of common variants
#			n.common.var <- length(intersect(encode.data[[i]], encode.data[[j]]))
			
			# get number of variants of reads i and j in overlaping region
#			n.var.overlap.i <- sum(encode.data[[i]]>=4*overlap.start & encode.data[[i]]<=4*overlap.end+3)
#			n.var.overlap.j <- sum(encode.data[[j]]>=4*overlap.start & encode.data[[j]]<=4*overlap.end+3)
			
				
 #                       if (n.var.overlap.i > min.overlap.var & n.var.overlap.j > min.overlap.var & 
#				n.common.var >= ceiling(n.var.overlap.i/2) & n.common.var >= ceiling(n.var.overlap.j/2)){
#                                read.pair.mat[i,j] <- TRUE
#                        }
#                }
#        }
#        read.pair.mat
	
#}


kmeans_reads <- function(encode.data, m5.data, centroid.seed, centroid.seed.range, min.cvg=20, min.overlap=200, max.iter=200)
{
	centroid.seed.update <- centroid.seed
	centroid.seed.range.update <- centroid.seed.range
	n.iter <- 0
	for (i in 1:max.iter){
		rl.core <- kmeans_reads_core(encode.data, m5.data, centroid.seed.update, centroid.seed.range.update,
					 min.cvg, min.overlap)
		if (identical(centroid.seed.update, rl.core$new.centroid))
			break
		centroid.seed.update <- rl.core$new.centroid
		centroid.seed.range.update <- rl.core$new.centroid.range
		n.iter <- n.iter + 1			
	}
	list(centroid = centroid.seed.update, centroid.range = centroid.seed.range.update, 
		centroid.read.id = rl.core$centroid.read.id, n.iter = n.iter)
						
}

kmeans_reads_core <- function(encode.data, m5.data, centroid.seed, centroid.seed.range, min.cvg=8, min.overlap=200, is.full.comp = FALSE)
{
	###--------- match reads to the centroid --------###
        if (length(encode.data)!=nrow(m5.data)){
                stop('encode.data and m5.data do not match.')
        }
	if (length(centroid.seed) != length(centroid.seed.range)) 
		stop('length(centroid.seed) != length(centroid.seed.range)')       

	group.id <- list() # 0 means consensus sequence (no variant)
	centroid.read.id <- lapply(centroid.seed, function(x) integer(0))
	for (i in 1:length(encode.data)){
                cur.jaccard <- rep(0, length(centroid.seed))
		cur.rate <- rep(0, length(centroid.seed))
		for (j in 1:length(centroid.seed)){
			overlap.start <- max(m5.data$tStart[i], centroid.seed.range[[j]][1])
			overlap.end <- min(m5.data$tEnd[i],  centroid.seed.range[[j]][2])	
			if (overlap.end - overlap.start + 1 < min.overlap)
				next
			
			cur.encode <- encode.data[[i]][encode.data[[i]] >= 4*overlap.start & encode.data[[i]] <= 4*overlap.end+3]
			cur.centroid <- centroid.seed[[j]][centroid.seed[[j]] >= 4*overlap.start & centroid.seed[[j]] <= 4*overlap.end+3]
			cur.common <- intersect(cur.encode, cur.centroid)
			cur.union <- union(cur.encode, cur.centroid)			

			# calculate jaccard index
			if (length(cur.union) > 0)
				cur.jaccard[j] <- length(cur.common) / length(cur.union)
			else
				cur.jaccard[j] <- 0
			
			# calculate proportion of shared variants in centroid
			if (length(cur.centroid) > 0)
				cur.rate[j] <- length(cur.common) / length(cur.centroid)
			else
				cur.rate[j] <- 0
				 
		}
		
		idx.on <- which(cur.rate>=0.5)
		if (length(idx.on) == 0){
                        group.id[[i]] <- 0
		}else{
			cur.jaccard.max <- max(cur.jaccard[idx.on])
			idx.jaccard.max <- idx.on[abs(cur.jaccard[idx.on] - cur.jaccard.max) <= 1e-12]
			group.id[[i]] <- idx.jaccard.max 
			for (j in 1:length(idx.jaccard.max))
				centroid.read.id[[idx.jaccard.max[j]]] <- c(centroid.read.id[[idx.jaccard.max[j]]], i)
		}
			
        }

	# update centroid accordingly
	new.centroid <- centroid.seed
	new.centroid.range <- centroid.seed.range
	for (i in 1:length(centroid.read.id)){
		if (length(centroid.read.id[[i]])==0){
			new.centroid[[i]] <- integer()
			new.centroid.range[[i]] <- c(0,0)
			break
		}
		cur.encode.data <- encode.data[centroid.read.id[[i]]]
		cur.m5.data <- m5.data[centroid.read.id[[i]], ]
		
		cur.max.encode <- max(4*m5.data$tEnd+3)
                cur.var.count <- integer(cur.max.encode + 1)
                cur.read.count <- integer(cur.max.encode + 1)

		# count number of variants and coverage for each locus		
		for (j in 1:length(cur.encode.data)){
			if (length(cur.encode.data[[j]]) == 0)
				next
			cur.var.count[cur.encode.data[[j]] + 1] <- cur.var.count[cur.encode.data[[j]] + 1] + 1
		}
		
		# count coverage for each locus
		for (j in 1:nrow(cur.m5.data)){
			cur.read.count[(4*cur.m5.data$tStart[j]+1):(4*cur.m5.data$tEnd[j]+3+1)] <- 
				cur.read.count[(4*cur.m5.data$tStart[j]+1):(4*cur.m5.data$tEnd[j]+3+1)] + 1
		}		
	
			
		# majority voting
		new.centroid[[i]] <- which(2*cur.var.count > cur.read.count) - 1

		# trim low coverage ends
		low.cvg <- which(cur.read.count < min.cvg) - 1
		new.centroid[[i]] <- setdiff(new.centroid[[i]], low.cvg)

		# remove variants outsides of new.centroid.range
		new.centroid[[i]] <- new.centroid[[i]][new.centroid[[i]] >= 4*new.centroid.range[[i]][1] & 
							new.centroid[[i]] <= 4*new.centroid.range[[i]][2]+3]
	
		# update range
		#high.cvg <- which(cur.read.count >= min.cvg) - 1
		#if (length(high.cvg) > 0){
		#	new.centroid.range[[i]] <- floor(range(high.cvg)/4)
		#}else{
		#	new.centroid.range[[i]] <- c(0,0)
		#}
		#new.centroid.range[[i]] <- floor(range(new.centroid[[i]])/4)
	}

	list(new.centroid=new.centroid, new.centroid.range=new.centroid.range, centroid.read.id=centroid.read.id)	
}




