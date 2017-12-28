#encode.data <- load.encodefile('./ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.encode.rdim.5000')
#m5.data <- load.m5file('./ERR752452_ERR690970_ERR1223274_ERR910547_ERR1588642.clean.m5.5000')
#k <- 200
#centroid.range <- c(m5.data$tStart[k], m5.data$tEnd[k])
#centroid <- encode.data[[k]]
#rl <- mat_fac_rank_1_core(encode.data, m5.data, centroid.range, centroid)
#rl <- mat_fac_rank_1(encode.data, m5.data, centroid.range, centroid, 10, 1000)
#rl <- mat_fac_rank_1(encode.data, m5.data, centroid.range, centroid, min.idx.on = 10, max.iter=100, is.full.comp = TRUE)
#rl <- mat_fac_rank_1(encode.data, m5.data, centroid.range, centroid, min.idx.on = 10, max.iter=100, is.full.comp = FALSE)


assemble_centroid <- function(centroid)
{
	centroid.start <- sapply(centroid, min)	
	centroid.sort <- centroid[order(centroid.start)]
	
}
trim_centroid <- function(encode.data, centroid, min.cvg = 10)
{
	centroid.trim <- centroid
	for (i in 1:length(centroid)){
		#print(i)
		if (length(centroid[[i]])==0)
			stop('length(centroid[[i]])==0')
		
		for (j in length(centroid[[i]]):1){
			n.var <- length(findReadsByLoci(encode.data, c(centroid[[i]][1],centroid[[i]][j])))				
			if (n.var >= min.cvg){
				centroid.trim[[i]] <- centroid[[i]][1:j]
				break
			}	
		}		
	}
	centroid.trim 
}

merge_centroid <- function(encode.data, m5.data, centroid, centroid.range.code, idx.on, min.rat=5)
{
	is.rm <- length(centroid)
	centroid.sim <- matrix(0, length(centroid), length(centroid))

	for (i in 1:length(centroid)){
		for (j in 1:length(centroid)){
			if (j == i)
				next
			if (min(centroid[[i]]) < min(centroid[[j]]) || max(centroid[[i]]) > max(centroid[[j]]) )
				next
		
			# only consider overlap region between centroids i an j
			#overlap.code.min <- max(centroid.range.code[[i]][1], centroid.range.code[[j]][1])
			#overlap.code.max <- min(centroid.range.code[[i]][2], centroid.range.code[[j]][2])
			
			overlap.code.min <- min(centroid[[i]])
			overlap.code.max <- max(centroid[[i]])
			if (overlap.code.min > overlap.code.max)
				next		
			
			centroid.i <- centroid[[i]][ centroid[[i]]>=overlap.code.min & centroid[[i]]<=overlap.code.max ]	
			centroid.j <- centroid[[j]][ centroid[[j]]>=overlap.code.min & centroid[[j]]<=overlap.code.max ]				
			idx.union <- sort(union(idx.on[[i]], idx.on[[j]])) 
			tStart.code <- 4*m5.data$tStart[idx.union]
			tEnd.code <- 4*m5.data$tEnd[idx.union] + 3

			# test different variant in centroid[i]
			diff.var.i <- setdiff(centroid.i, centroid.j)
			n.var.i <- numeric(length(diff.var.i))
			n.reads.i <- numeric(length(diff.var.i))
			if (length(n.var.i)>0){
				for (k in 1:length(n.var.i)){
					n.var.i[k] <- length(findReadsByLoci(encode.data[idx.union], diff.var.i[k]))
					n.reads.i[k] <- sum(tStart.code<=min(diff.var.i[k]) & tEnd.code>=max(diff.var.i[k]))
				}
			}
	
			# test different variant in centroid[j]
			diff.var.j <- setdiff(centroid.j, centroid.i)
			n.var.j <- numeric(length(diff.var.j))
                        n.reads.j <- numeric(length(diff.var.j))
			if (length(n.var.j)>0){
                        	for (k in 1:length(n.var.j)){
                                	n.var.j[k] <- length(findReadsByLoci(encode.data[idx.union], diff.var.j[k]))
                                	n.reads.j[k] <- sum(tStart.code<=min(diff.var.j[k]) & tEnd.code>=max(diff.var.j[k]))
                        	}
			}
			
		}
	}		
}

mat_fac_old <- function(encode.data, m5.data, min.cvg = 100,  min.idx.on = 10, max.iter=1)
{
	# variants set 
	var.set <- sort(unique(unlist(encode.data)))

	# scan var.set and find largest windows covered by # of reads >= min.cvg
	is.over = FALSE
	i <- 1
	#range.start <- var.set[i]
	#range.end <- var.set[i]
	while(TRUE){
		for (j in i:length(var.set)){
			if (j+1 > length(var.set)){
				is.over <- TRUE
				break
			}
			n.read.cover <- sum(m5.data$tStart<=var.set[i]  & m5.data$tEnd>=var.set[j+1] )
			if (n.read.cover < min.cvg)
				break	
		}
		idx <- which(m5.data$tStart<=var.set[i] & m5.data$tEnd>=var.set[j])
		centroid.idx <- idx[which.max( (m5.data$tEnd-m5.data$tStart)[idx] )]

		centroid <- encode.data[[centroid.idx]]
		centroid.range <- c(floor(var.set[i]/4), floor(var.set[j]/4))
	}
}

mat_fac_rm_redudant <- function(centroid)
{
	is.redundant <- rep(FALSE, length(centroid))
	for (i in 1:length(centroid)){
        	for (j in 1:length(centroid)){
                	if (j==i)
                        	next
                	overlap.start <- min(centroid[[i]])
                	overlap.end <- max(centroid[[i]])

                	centroid.i <- centroid[[i]]
                	centroid.j <- centroid[[j]][centroid[[j]]>=overlap.start & centroid[[j]]<=overlap.end]
                	if (identical(centroid.i, centroid.j) & length(centroid[[i]]) <= length(centroid[[j]]) ){
                        	is.redundant[i] <- TRUE
                        	break
                	}
        	}
	}
	centroid[!is.redundant]
}

cal.dist <- function(encode.1, range.1, encode.2, range.2)
{
	range.overlap <- c(max(range.1[1], range.2[1]), min(range.1[2], range.2[2]))
	if (range.overlap[1] > range.overlap[2])
		return(NaN)
	encode.1.overlap <- encode.1[encode.1>=4*range.overlap[1] & encode.1<=4*range.overlap[2]+3]
	encode.2.overlap <- encode.2[encode.2>=4*range.overlap[1] & encode.2<=4*range.overlap[2]+3]		
	
	encode.dist <- length(setdiff(encode.1.overlap, encode.2.overlap)) + length(setdiff(encode.2.overlap, encode.1.overlap))
	encode.dist <- encode.dist / (range.overlap[2] - range.overlap[1])	
	encode.dist	
}

# use uncontained reads as seed centroids
mat_fac <- function(encode.data, m5.data, centroid.seed, centroid.seed.range, min.cvg = 20, min.idx.on = 10, max.iter=100, is.full.comp = FALSE, is.overhang=TRUE)
{
	# check if encode.data matches m5.data
        if (length(encode.data) != nrow(m5.data))
                stop('length(encode.data) != nrow(m5.data)')
	# check if centroid.seed and centroid.seed.range compatible
	if (length(centroid.seed) != length(centroid.seed.range))
		stop('length(centroid.seed) != length(centroid.seed.range)')

	# conduct rank 1 decomposition for each centroid.seed
	mat.fac.rl <- list()
	for (i in 1:length(centroid.seed)){
		cat('centroid.seed', i, '\n')
		# rank 1 matrix facterization 
		if (length(centroid.seed[[i]])==0)
			next

		cur.rl <- mat_fac_rank_1(encode.data, m5.data, centroid.seed.range[[i]], centroid.seed[[i]], min.cvg, min.idx.on,
                                                max.iter, is.full.comp, is.overhang)
		if (is.null(cur.rl))
			next
		if (length(cur.rl$idx.on) >= min.idx.on){
			mat.fac.rl[[i]] <- cur.rl 
		}

	}

	mat.fac.rl			
}
# use the read with largest minial distance to the existing centroids as the new seed centroid
mat_fac_mindist <- function(encode.data, m5.data, min.cvg = 20, min.idx.on = 10, max.iter=100, is.full.comp = FALSE, is.overhang=TRUE)
{
	# check if encode.data matches m5.data
        if (length(encode.data) != nrow(m5.data))
                stop('length(encode.data) != nrow(m5.data)')
	
	# set origin to be the first centroid and select the read that have largest minial distance to all existing centroid
	# to be the next centroid
	centroid.list <- list(integer(0))
	centroid.range.list <- list(c(min(m5.data$tStart), max(m5.data$tEnd)))
	idx.on.all <- integer(0)
        idx.off.all <- 1:length(encode.data)
        mat.fac.rl <- list()
        k <- 1
	while(TRUE){
		read.dist <- rep(1, length(idx.off.all))
		# scan each read
		for (i in 1:length(idx.off.all)){
			cur.encode <- encode.data[[ idx.off.all[i] ]]
			cur.range <- c(m5.data$tStart[idx.off.all[i]], m5.data$tEnd[idx.off.all[i]])
			# find mimial distance to the existing centroids
			for (j in 1:length(centroid.list)){
				cur.dist <- cal.dist(cur.encode, cur.range, centroid.list[[j]], centroid.range.list[[j]])	
				if (cur.dist < read.dist[i])
					read.dist[i] <- cur.dist
			}
		}
		idx.max <- idx.off.all[which.max(read.dist)]
		
		# rank 1 matrix facterization 
		cur.centroid <- encode.data[[idx.max]]
		cur.centroid.range <- c(m5.data$tStart[idx.max], m5.data$tEnd[idx.max])
	
		if (length(cur.centroid)==0){
			idx.on.all <- union(idx.on.all, idx.max)
		}else{
			cur.rl <- mat_fac_rank_1(encode.data, m5.data, cur.centroid.range, cur.centroid, min.cvg, min.idx.on,
                                                max.iter, is.full.comp, is.overhang)	
			
			# if size of idx.on of current decomposition is larger than min.idx.on
                	# trim centroid, refacterize and update idx.on.all and idx.off.all
                	if (length(cur.rl$idx.on) >= min.idx.on){
                		#cur.centroid.trim <- trim_centroid(encode.data, list(cur.rl$new.centroid), min.cvg)[[1]]
                        	#cur.centroid.trim.range <- c(floor(min(cur.centroid.trim)/4), floor(max(cur.centroid.trim)/4))
                        	#cur.rl <- mat_fac_rank_1(encode.data, m5.data, cur.centroid.trim.range, cur.centroid.trim,
                                #			min.idx.on, max.iter, is.full.comp)
                        	if (all(cur.rl$idx.on %in% idx.on.all)){
                        		idx.on.all <- union(idx.on.all, idx.max)
					#next
				}

                        	mat.fac.rl[[k]] <- cur.rl
                        	k <- k + 1
               		}
		}
		idx.on.all <- union(idx.on.all, cur.rl$idx.on)
                idx.off.all <- setdiff(1:length(encode.data), idx.on.all)
		cat("idx.off.all size:", length(idx.off.all), '\n')
		if (length(idx.off.all)==0)
			break		
	}
	mat.fac.rl		
}

mat_fac_long_reads_first <- function(encode.data, m5.data, min.cvg = 10, min.idx.on = 10, max.iter=100, is.full.comp = TRUE)
{	
	# check if encode.data matches m5.data
	if (length(encode.data) != nrow(m5.data))
		stop('length(encode.data) != nrow(m5.data)')
	
	# sort reads according to mapped length
	idx.order <- order(m5.data$tEnd - m5.data$tStart, decreasing=TRUE)		
	encode.data.sort <- encode.data[idx.order]	
	m5.data.sort <- m5.data[idx.order, ]
	
	# use the longest read that is not in idx.on.all
	idx.on.all <- integer(0)
	idx.off.all <- 1:length(encode.data.sort)
	mat.fac.rl <- list()
	k <- 1	
	while(TRUE){
		cat("idx.off.all size:", length(idx.off.all), '\n')
		if (length(idx.off.all)==0)
			break
		for (i in 1:length(idx.off.all)){
			cur.centroid <- encode.data.sort[[ idx.off.all[i] ]]
			cur.centroid.range <- c(m5.data.sort$tStart[ idx.off.all[i] ], m5.data.sort$tEnd[ idx.off.all[i] ])
			if (length(cur.centroid) == 0)
				next
			
				
			# rank 1 decomposition 
			cur.rl <- mat_fac_rank_1(encode.data.sort, m5.data.sort, cur.centroid.range, cur.centroid, min.idx.on, 
						max.iter, is.full.comp)
			
			# if size of idx.on of current decomposition is larger than min.idx.on 
			# trim centroid, refacterize and update idx.on.all and idx.off.all
			if (length(cur.rl$idx.on) >= min.idx.on){
				cur.centroid.trim <- trim_centroid(encode.data.sort, list(cur.rl$new.centroid), min.cvg)[[1]]
				cur.centroid.trim.range <- c(floor(min(cur.centroid.trim)/4), floor(max(cur.centroid.trim)/4))
				cur.rl <- mat_fac_rank_1(encode.data.sort, m5.data.sort, 
							cur.centroid.trim.range, cur.centroid.trim, 
							min.idx.on, max.iter)
				if (all(cur.rl$idx.on %in% idx.on.all))
					next
				
				mat.fac.rl[[k]]	<- cur.rl
				k <- k + 1
				idx.on.all <- union(idx.on.all, cur.rl$idx.on)
				idx.off.all <- setdiff(1:length(encode.data.sort), idx.on.all)
				break
			}
		}
		cat('i =', i, '\n')
		if (i == length(idx.off.all))
			break
	}
	mat.fac.rl	
}

mat_fac_rank_1 <- function(encode.data, m5.data, centroid.range, centroid, min.cvg =20,min.idx.on = 10, max.iter=100,
			 is.full.comp = FALSE, is.overhang=TRUE)
{
	old.centroid <- centroid
	n.iter <- 0
	for (i in 1:max.iter){
		#print(i)
		rl <- mat_fac_rank_1_core(encode.data, m5.data, centroid.range, old.centroid, min.cvg, is.full.comp, is.overhang)
		if (is.null(rl)){
			return(NULL)
		}
		n.iter <- n.iter + 1
		if (length(rl$idx.on) < min.idx.on)
			break;
			
		if (identical(as.integer(rl$new.centroid), as.integer(old.centroid)))
                        break
                
                old.centroid <- rl$new.centroid
                old.idx.on <- rl$idx.on
		old.idx.off <- rl$idx.off
		centroid.range <- rl$new.centroid.range

	}
	rl$n.iter <- n.iter
	rl
}

mat_fac_rank_1_core <- function(encode.data, m5.data, centroid.range, centroid, min.cvg = 20, is.full.comp = FALSE, is.overhang = TRUE)
{
	###--------- match reads to the centroid --------###
	if (length(encode.data)!=nrow(m5.data)){
		stop('encode.data and m5.data do not match.')
	}
	is.on <- rep(FALSE, length(encode.data))
	for (i in 1:length(encode.data)){
		if (is.full.comp){
			common.var <- intersect(encode.data[[i]], centroid)
			if (length(common.var) >= ceiling(length(centroid)/2 ))
				is.on[i] <- TRUE
		}else{
			common.var <- intersect(encode.data[[i]], centroid)
			n.centroid.overlap <- sum(centroid>=4*m5.data$tStart[i] & centroid<=4*m5.data$tEnd[i]+3)
			if (length(common.var)>n.centroid.overlap)
				stop('length(common.var)>n.centroid.overlap')
			if (length(common.var) > ceiling(n.centroid.overlap/2 ) & 
				n.centroid.overlap > 0)
                                is.on[i] <- TRUE
		}
	}
	idx.on <- which(is.on)
	idx.off <- which(!is.on)
	
	###--------- update centroid -------###
	# use the binary matrix to represent reads
	centroid.range.code <- centroid.range
	
	if (is.overhang){
		centroid.range[1] <- min(m5.data$tStart[idx.on])
		centroid.range[2] <- max(m5.data$tEnd[idx.on])		
	}
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
	
	# trim new.centroid 
	cvg.profile <- apply(encode.data.mat.reads, 2, sum)
	idx.true <- which(cvg.profile >= min.cvg)
	if (length(idx.true)==0){
		return(NULL)
	}
	new.centroid.range <- floor((range(idx.true) + centroid.range.code[1] - 1) / 4)
	new.centroid <- new.centroid[new.centroid>=4*new.centroid.range[1] & new.centroid<=4*new.centroid.range[2]+3]
	list(idx.on=idx.on, idx.off=idx.off, new.centroid=new.centroid, new.centroid.range=new.centroid.range)
}



