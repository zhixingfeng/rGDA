library(dequer)
check.contain <- function(rl.cluster){
	contain.mat <- matrix(0, length(rl.cluster), length(rl.cluster))	
	for (i in 1:length(rl.cluster)){
		for (j in 1:length(rl.cluster)){
			if (i==j) next
			overlap.start <- max(rl.cluster[[i]]$centroid.range[1], rl.cluster[[j]]$centroid.range[1])
			overlap.end <- min(rl.cluster[[i]]$centroid.range[2], rl.cluster[[j]]$centroid.range[2])
			overlap.len <- overlap.end - overlap.start
			if (rl.cluster[[i]]$centroid.range[1] < rl.cluster[[j]]$centroid.range[1] | 
				rl.cluster[[i]]$centroid.range[2] > rl.cluster[[j]]$centroid.range[2])	
				next
			centroid.i <- rl.cluster[[i]]$centroid
			centroid.j <- rl.cluster[[j]]$centroid
			centroid.j <- centroid.j[centroid.j >= 4*overlap.start & centroid.j <= 4*overlap.end+3]
			
			if (identical(centroid.i, centroid.j))
				contain.mat[i,j] <- 1
		}
	}
	is.contained <- apply(contain.mat, 1, function(x) any(x==1))
	is.contained
}


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
			max.iter = 50, is.trim=TRUE, is.full.comp = TRUE)
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
                                        min.overlap = min.overlap, is.full.comp = is.full.comp)
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

assemble.core <- function(encode.data, m5.data, centroid, centroid.range, min.overlap=200, is.full.comp = TRUE)
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
	adj.mat <- matrix(0, length(rl.cluster), length(rl.cluster))
	for (i in 1:length(rl.cluster)){
		for (j in 1:length(rl.cluster)){
			if (i==j) next
			if (rl.cluster[[i]]$centroid.range[2] > rl.cluster[[j]]$centroid.range[2])
				next
			overlap.start <- max(rl.cluster[[i]]$centroid.range[1], rl.cluster[[j]]$centroid.range[1])
			overlap.end <- min(rl.cluster[[i]]$centroid.range[2], rl.cluster[[j]]$centroid.range[2])
			
			if (overlap.start > overlap.end) next			
			centroid.i <- rl.cluster[[i]]$centroid[rl.cluster[[i]]$centroid >= 4*overlap.start &
							rl.cluster[[i]]$centroid <= 4*overlap.end+3]
			centroid.j <- rl.cluster[[j]]$centroid[rl.cluster[[j]]$centroid >= 4*overlap.start &
                                                        rl.cluster[[j]]$centroid <= 4*overlap.end+3]
			if (identical(centroid.i, centroid.j) & length(centroid.i)>=1)
				adj.mat[i,j] <- 1
		}
	}
	adj.mat
}


link.scafold.bak <- function(rl.cluster)
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
				if (identical(cur.haplotype, cur.centroid) & length(cur.haplotype) >= 1){
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

merge.scafold <- function(adj.mat)
{
	start.nodes <- which(apply(adj.mat, 2, sum) == 0)
        adj.list <- apply(adj.mat, 1, function(x) which(x==1))
	path.list <- list()
	while(TRUE){
		if (length(start.nodes) == 0 | is.null(start.nodes)){
                        break
		}else{
			rl.core <- merge.scafold.core(adj.list, start.nodes)
			path.list <- c(path.list, rl.core$path.list)
			start.nodes <- rl.core$new.start.nodes
		}
	}
	path.list
}

merge.scafold.core <- function(adj.list, start.nodes)
{
	path.list <- list()
	new.start.nodes <- list()
	n.path <- 0
	for (cur.start.node in start.nodes){
		n.path <- n.path + 1
		cur.node <- cur.start.node
		cur.path <- stack()
		while(TRUE){
			push(cur.path, cur.node)
			# if number of daughter nodes is 1, then get into deeper
			if (length(adj.list[[cur.node]]) == 1){
				cur.node <- adj.list[[cur.node]][[1]]
				next
			}
			
			# if no daughter nodes, stop
			if (length(adj.list[[cur.node]]) == 0){
				new.start.nodes[[n.path]] <- integer(0)
				path.list[[n.path]] <- rev(unlist(as.list(cur.path)))
                                break
                        }
			
			# if more than 1 daughter nodes, stop and add the daughter nodes to new start.nodes
			if (length(adj.list[[cur.node]]) > 1){
				new.start.nodes[[n.path]] <- adj.list[[cur.node]]
				path.list[[n.path]] <- rev(unlist(as.list(cur.path)))
                                break
                        }
		}
	}
		
	list(path.list = path.list, new.start.nodes = sort(unique(unlist(new.start.nodes))))
}

merge.scafold.extend <- function(adj.mat)
{
	start.nodes <- which(apply(adj.mat, 2, sum) == 0)
        adj.list <- apply(adj.mat, 1, function(x) which(x==1))
        path.list <- list()
	
	while (length(start.nodes) > 0){
		cur.rl <- merge.scafold.extend.core(adj.list, start.nodes)
		start.nodes <- cur.rl$new.start.nodes
		path.list <- c(path.list, cur.rl$path.list)
	}
	path.list
}

merge.scafold.extend.core <- function(adj.list, start.nodes)
{
	path.list <- list()
	new.start.nodes <- integer()
        for (cur.start.node in start.nodes){
		# follow unambigous apth from cur.node
        	rl <- merge.scafold.core(adj.list, cur.start.node)
		
		# number of new start nodes should never be 1
		if (length(rl$new.start.nodes) == 1)
			stop('length(rl$new.start.nodes) == 1')
		# if hit the end, stop
		if (length(rl$new.start.nodes) == 0){
			path.list <- c(path.list, rl$path.list)
			next
		}

		# if hit divided paths, follow each of them
		if (length(rl$new.start.nodes) > 1){
			rl.new <- merge.scafold.core(adj.list, rl$new.start.node)
			for (i in 1:length(rl.new$path.list)){
				rl.new$path.list[[i]] <- c(rl$path.list[[1]], rl.new$path.list[[i]])
			}
			path.list <- c(path.list, rl.new$path.list)
			
			add.new.start.nodes <- rl.new$new.start.nodes[sapply(adj.list[rl.new$new.start.nodes],length) > 0]
			new.start.nodes <- c(new.start.nodes, add.new.start.nodes)
		}
	}
	
	list(path.list = path.list, new.start.nodes = sort(unique(new.start.nodes)))
}



