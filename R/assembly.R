assemble.trim.centroid <- function(m5.data.sub, centroid, centroid.range, min.cvg = 20)
{
	# if centroid has no variants quit
	if (length(centroid)==0)
		return(NULL)

	# if centroid has only 1 variant, trim it directly
	rl <- list()
	if (length(centroid)==1){
		rl[[1]] <- assemble.trim.centroid.core(m5.data.sub, centroid, centroid.range, min.cvg)
		return(rl)
	}
	
	# if centroid has multiple variants, scan it and trim. 
	cut.points.left <- floor(centroid[-length(centroid)]/4)
	cut.points.right <- floor(centroid[-1]/4)
	cut.points <- floor((cut.points.right + cut.points.left)/2) 
	cut.points <- c(centroid.range[1], cut.points)
	for (i in 1:length(cut.points)){
		cur.centroid.range <- c(cut.points[i], centroid.range[2])
		cur.centroid <- centroid[centroid >= 4*cur.centroid.range[1] & centroid <= (4*cur.centroid.range[2]+3)]
		cur.rl <- assemble.trim.centroid.core(m5.data.sub, cur.centroid, cur.centroid.range, min.cvg)
		if (is.null(cur.rl))
			break
		rl[[i]] <- assemble.trim.centroid.core(m5.data.sub, cur.centroid, cur.centroid.range, min.cvg)
		if (rl[[i]]$centroid.range.trim[2] == centroid.range[2])
			break
	}
	rl
}
assemble.trim.centroid.core <- function(m5.data.sub, centroid, centroid.range, min.cvg = 20)
{
	cvg.profile <- rep(0, 4*max(m5.data.sub$tEnd) + 3 + 1)	
	idx <- which(m5.data.sub$tStart <= centroid.range[1] & m5.data.sub$tEnd > centroid.range[1])
	for (i in idx){
		cur.loci <- (4*m5.data.sub$tStart[i] + 1) : (4*m5.data.sub$tEnd[i] + 3 + 1)
		cvg.profile[cur.loci] <- cvg.profile[cur.loci] + 1
	}
	
	if (sum(cvg.profile >= min.cvg) == 0)
		return(NULL)

	centroid.range.upper <- floor((max(which(cvg.profile >= min.cvg)) - 1) / 4)
	if (centroid.range[1] >= centroid.range.upper)
		return(NULL)
	centroid.range.trim <- c(centroid.range[1], min(centroid.range.upper, centroid.range[2]))
        centroid.trim <- centroid[centroid >= 4*centroid.range.trim[1] & centroid <= 4*centroid.range.trim[2] + 3]
        if (length(centroid.trim)==0)
		return(NULL)
	list(centroid.trim = centroid.trim, centroid.range.trim = centroid.range.trim)
}

assemble <- function(encode.data, m5.data, centroid, centroid.range, min.idx.on = 20, min.cvg = 20, min.overlap = 200, 
		max.iter = 50, is.trim=TRUE)
{
	#------------- check input data ------------#
        if (length(encode.data)!=nrow(m5.data)){
                stop('encode.data and m5.data do not match.')
        }

	if (length(centroid)==0)
		return(NULL)
	if (any(centroid < 4*centroid.range[1] | centroid > (4*centroid.range[2]+3)))
		stop('centroid and centroi.range do not match')	

	#------------- update centroid -----------#
	encode.data.sub <- encode.data
        m5.data.sub <- m5.data
        centroid.update <- centroid
        centroid.range.update <- centroid.range
	for (i in 1:max.iter){
		rl.core <- assemble.core(encode.data.sub, m5.data.sub,
                                        centroid.update, centroid.range.update,
                                        min.overlap = min.overlap)
		if (length(rl.core$centroid)==0 | identical(rl.core$centroid, centroid.update) | length(rl.core$idx.on) < min.idx.on)	
			break
		centroid.update <- rl.core$centroid
		centroid.range.update <- rl.core$centroid.range
	}
	rl.core[['n.iter']] <- i

	#------------- trim centroid -------------#
	if (is.trim){	
		if (length(rl.core$idx.on) < min.cvg)
			return(rl.core)
		rl.trim <- assemble.trim.centroid(m5.data[rl.core$idx.on,], rl.core$centroid, rl.core$centroid.range, min.cvg)
		rl.core[['trim']] <- rl.trim
	}
	rl.core[['centroid.init']] <- centroid
	rl.core
}

assemble.left <- function(encode.data, m5.data, centroid, centroid.range, min.cvg = 20, min.overlap = 200, max.iter = 50)
{
	 #------------- check input data ------------#
        if (length(encode.data)!=nrow(m5.data)){
                stop('encode.data and m5.data do not match.')
        }

	idx.on <- 1:length(encode.data)
	for (i in 1:max.iter){
		# trim centroid
		idx <- idx.on[m5.data$tStart[idx.on] <= centroid.range[1] & m5.data$tEnd[idx.on] > centroid.range[1]]
		
		cvg.profile <- rep(0, 4*centroid.range[2] + 3 - 4*centroid.range[1] + 1)
		for (j in idx){
			cur.start <- max(m5.data$tStart[j], centroid.range[1])
			cur.end <- min(m5.data$tEnd[j], centroid.range[2])
			cur.loci <- (4*(cur.start - centroid.range[1]) + 1) : (4*(cur.end - centroid.range[1]) + 3 + 1)
			cvg.profile[cur.loci] <- cvg.profile[cur.loci] + 1
		}
		if (sum(cvg.profile >= min.cvg) == 0)
			break
		centroid.range.upper <- floor((max(which(cvg.profile >= min.cvg)) - 1 + 4*centroid.range[1]) / 4)
		centroid.range <- c(centroid.range[1], min(centroid.range.upper, centroid.range[2]))
		centroid <- centroid[centroid >= 4*centroid.range[1] & centroid <= 4*centroid.range[2] + 3]	
		if (length(centroid)==0)
			break

		# update centroid
		rl.core <- assemble.core(encode.data, m5.data, centroid, centroid.range, min.overlap = min.overlap)
		idx.on <- rl.core$idx.on
		if (identical(centroid, rl.core$centroid))
			break
		centroid <- rl.core$centroid
	}

	list(centroid = centroid, centroid.range = centroid.range, idx.on = idx.on, n.iter = i)
}

assemble.core <- function(encode.data, m5.data, centroid, centroid.range, min.overlap=200)
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
                n.centroid.overlap <- sum(centroid>=4*m5.data$tStart[i] & centroid<=4*m5.data$tEnd[i]+3)
                if (length(common.var)>n.centroid.overlap)
                	stop('length(common.var)>n.centroid.overlap')
                if (length(common.var) >= ceiling(n.centroid.overlap/2) & n.centroid.overlap > 0)
                        is.on[i] <- TRUE
        }
        idx.on <- which(is.on)	

	#------------- update centroid -------------#
	new.centroid <- centroid
	cur.encode.data <- encode.data[idx.on]
	cur.m5.data <- m5.data[idx.on,]
	
	cur.max.encode <- max(4*m5.data$tEnd+3)
        cur.var.count <- integer(cur.max.encode + 1)
        cur.read.count <- integer(cur.max.encode + 1)

	# count number of variants and coverage for each locus          
        for (i in 1:length(cur.encode.data)){
		if (length(cur.encode.data[[i]]) == 0)
			next
                cur.var.count[cur.encode.data[[i]] + 1] <- cur.var.count[cur.encode.data[[i]] + 1] + 1
        }

        # count coverage for each locus
        for (i in 1:nrow(cur.m5.data)){
                cur.read.count[(4*cur.m5.data$tStart[i]+1):(4*cur.m5.data$tEnd[i]+3+1)] <-
                	cur.read.count[(4*cur.m5.data$tStart[i]+1):(4*cur.m5.data$tEnd[i]+3+1)] + 1
        }

        # majority voting
        new.centroid <- which(cur.var.count>=ceiling(cur.read.count/2) & cur.read.count>0) - 1
	
	# restrict new.centroid to centroid.range
	new.centroid <- new.centroid[new.centroid >= 4*centroid.range[1] & new.centroid <= 4*centroid.range[2]+3]

	list(centroid = new.centroid, centroid.range = centroid.range,idx.on = idx.on)
}


link.scafold <- function(rl.cluster)
{
	haplotypes <- list()
	haplotypes.range <- list()
	for (i in 1:length(rl.cluster)){
		print(i)
		if (length(haplotypes)==0){
			haplotypes[[1]] <- rl.cluster[[i]]$centroid
			haplotypes.range[[1]] <- rl.cluster[[i]]$centroid.range
			next
		}
		
		is.new <- TRUE
		for (j in 1:length(haplotypes)){
			overlap.start <- max(haplotypes.range[[j]][1], rl.cluster[[i]]$centroid.range[1])
			overlap.end <- min(haplotypes.range[[j]][2], rl.cluster[[i]]$centroid.range[2])
			if (overlap.start < overlap.end ){
				cur.haplotype <- haplotypes[[j]][haplotypes[[j]] >= 4*overlap.start & 
								haplotypes[[j]] <= 4*overlap.end+3]
				cur.centroid <- rl.cluster[[i]]$centroid[rl.cluster[[i]]$centroid >= 4*overlap.start &
								rl.cluster[[i]]$centroid <= 4*overlap.end+3]
				if (identical(cur.haplotype, cur.centroid)){
					haplotypes[[j]] <- sort(unique(c(haplotypes[[j]], rl.cluster[[i]]$centroid)))
					haplotypes.range[[j]][1] <- min(haplotypes.range[[j]][1], rl.cluster[[i]]$centroid.range[1])
					haplotypes.range[[j]][2] <- max(haplotypes.range[[j]][2], rl.cluster[[i]]$centroid.range[2])
					is.new <- FALSE
				}
			}

		}
		if (is.new){
			haplotypes <- c(haplotypes, list(rl.cluster[[i]]$centroid))
			haplotypes.range <- c(haplotypes.range, list(rl.cluster[[i]]$centroid.range))
		}
	}
	list(haplotypes=haplotypes, haplotypes.range=haplotypes.range)
}


link.scafold.bak <- function(rl.cluster)
{
	node.all <- lapply(rl.cluster, function(x) x$centroid)
	node.all <- sort(unique(unlist(node.all)))
	adj.mat <- matrix(0, length(node.all), length(node.all))
	rownames(adj.mat) <- node.all
	colnames(adj.mat) <- node.all
	for (i in 1:length(rl.cluster)){
		print(i)
		for (j in 1:length(rl.cluster[[i]]$trim)){
			idx <- match(rl.cluster[[i]]$trim[[j]]$centroid.trim, node.all)
			if (length(idx) < 2)
				next
				#stop('length(idx) < 2')
			for (k in 1:(length(idx)-1)){
				adj.mat[idx[k], idx[k+1]] <- 1
			}
		}
	}
	adj.mat	
}


