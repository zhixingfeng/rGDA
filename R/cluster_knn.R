cluster.knn <- function(encode.data, m5.data, topn = 20, step.size = 20, max.nn = 200)
{
	cat('calculate pairwise distance\n')
	encode.mat <- conv_encode_to_mat(encode.data) + 0
	reads.dist <- as.matrix(dist(encode.mat))^2

	cat('clustering :\n')
	neighbor.id.list <- list()
	cons.list <- list()	
	for (i in 1:nrow(reads.dist)){
		cat('read :',i,'\n')
		cur.order <- order(reads.dist[i,])
		n.size <- topn
		neighbor.id <- NULL
		cur.neighbor.id <- NULL
		while(n.size <= max.nn){
			cur.neighbor.id <- cur.order[1:n.size]
			cur.cons <- get_consensus(encode.data[cur.neighbor.id], m5.data[cur.neighbor.id,], rm.del = FALSE)
			if (any(cur.cons$prop>=0.2 & cur.cons$prop<=0.7,na.rm = TRUE)){
				if (n.size > topn){
					neighbor.id <- cur.order[1:(n.size - step.size)]
					cur.cons <- get_consensus(encode.data[neighbor.id], m5.data[neighbor.id,], rm.del = FALSE)
				}
				break
			}
			n.size <- n.size + step.size
		}

		if (n.size > max.nn)
			neighbor.id <- cur.neighbor.id 

		neighbor.id.list[[i]] <- neighbor.id
		cons.list[[i]] <- cur.cons
	}
	
	list(neighbor.id.list = neighbor.id.list, cons.list = cons.list)	

}


