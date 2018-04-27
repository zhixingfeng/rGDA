slide_win_clustering_core <- function(encode.data.win, m5.data.win, tmp.dir, min.jaccard = 0.8, prefix = 'win_xx_xx')
{
	# save encode.data and m5.data
	encode.file.win <- paste(tmp.dir, '/', prefix, '.encode', sep = '')
	m5.file.win <- paste(tmp.dir, '/', prefix, '.m5', sep = '')

	save.encodefile(encode.data.win, encode.file.win)
	save.m5file(m5.data.win, m5.file.win)
	
	# calculate jaccard index
	jaccard.file <- paste(tmp.dir, '/', prefix, '.jaccard', sep = '')
	cmd <- paste('igda jaccard -m 0.8', encode.file.win, m5.file.win, jaccard.file)
	print(cmd); system(cmd)

	# load jaccard index and build graph
	jaccard.index <- read.table(jaccard.file, header=FALSE, sep=',')
	g.dat <- jaccard.index[,1:3]
	colnames(g.dat)[3] <- 'weight'
	jaccard.index.graph <- graph.data.frame(g.dat, directed = FALSE)
	
	# infomap clustering
	cl.infomap <- cluster_infomap(jaccard.index.graph)
	
	# consensus
	cl.infomap.gp <- split(as.integer(cl.infomap$names)+1, cl.infomap$membership)
	rl.cons <- list()
	for (i in 1:length(cl.infomap.gp)){
		rl.cons[[i]] <- get_consensus(encode.data.win[cl.infomap.gp[[i]]],  m5.data.win[cl.infomap.gp[[i]],])
	}
	rl.cons
}

get_consensus <- function(encode.data.gp, m5.data.gp, rm.del = TRUE)
{
	max.encode <- max(4*m5.data.gp$tEnd+3)
        var.count <- rep(0,max.encode + 1)
        read.count <- rep(0,max.encode + 1)			

	# count number of variants and coverage for each locus
        for (i in 1:length(encode.data.gp)){
                if (length(encode.data.gp[[i]]) == 0)
                        next
                var.count[encode.data.gp[[i]] + 1] <- var.count[encode.data.gp[[i]] + 1] + 1
        }

        # count coverage for each locus
	if (rm.del){
		read.count <- pileup_reads_count(m5.data.gp)
	}else{
	        for (i in 1:nrow(m5.data.gp)){
                	read.count[(4*m5.data.gp$tStart[i]+1):(4*m5.data.gp$tEnd[i]+3+1)] <-
                	        read.count[(4*m5.data.gp$tStart[i]+1):(4*m5.data.gp$tEnd[i]+3+1)] + 1
        	}
	}

	# majority vote 
	cons.seq <- which(var.count>=ceiling(read.count/2) & read.count>0) - 1

	list(cons.seq = cons.seq, var.count = var.count, read.count = read.count, prop = var.count / read.count)
}

get_consensus_recode <- function(encode.data.gp, m5.data.gp, min.read.count = 5)
{
	max.encode <- max(4*m5.data.gp$tEnd+3)
	pu.var <- pileup_var(encode.data.gp, max.encode)
	pu.read <- pileup_reads(m5.data.gp)
	
	var.count <- sapply(pu.var, length)
	read.count <- sapply(pu.read, length)
	
	prop <- mapply(function(x,y) length(intersect(x,y)) / length(y), x = pu.var, y = pu.read)
	#prop <- mapply(function(x,y) length(x) / length(union(x,y)), x = pu.var, y = pu.read)

	cons.seq <- which(prop >0.5 & read.count >= min.read.count) - 1

	list(cons.seq = cons.seq, var.count = var.count, read.count = read.count, prop = prop)
}

get_cons_bymovie <- function(encode.data, m5.data, min.read.count = 5)
{
        if (!all(validate.data(encode.data, m5.data)))
                stop('!all(validate.data(encode.data, m5.data))')
        if (!all(check.m5(m5.data)))
                stop('!all(check.m5(m5.data))')
        
        movies <- get.movie(m5.data)
        encode.data.movie <- split(encode.data, movies)
        m5.data.movie <- split(m5.data, movies)

        cons.list <- list()
        for (i in 1:length(encode.data.movie))
                cons.list[[i]] <- get_consensus_recode(encode.data.movie[[i]], m5.data.movie[[i]], min.read.count)
        cons.list
}










