### abandon
#library(dequer)

unconvolve.scafold <- function(adj.mat)
{
	is.start <- apply(adj.mat, 2, sum) == 0
	is.end <- apply(adj.mat, 1, sum) == 0
	
	adj.list <- apply(adj.mat, 1, function(x) which(x==1))
	adj.list.visited <- lapply(adj.list, function(x) rep(FALSE, length(x)))
	node.next <- rep(1, length(adj.list))
	node.prev <- rep(0, length(adj.list))
	path.list <- list()
	n.path <- 0
	for (start.node in which(is.start)){
		# deep first search 
		is.confused <- FALSE
		cur.node <- start.node
		cur.path <- stack()
		push(cur.path, cur.node)
		
		while (TRUE){	
			#print(cur.node); readline()
			# the start.node has visited and trace back to its previous node, which is 0
                        if (cur.node == 0)
                                break
			# trace back if a node has no daughter node or all its daughter nodes are visited
			is.end <- length(adj.list.visited[[cur.node]]) == 0
			is.visited <- length(adj.list.visited[[cur.node]]) > 0 & all(adj.list.visited[[cur.node]]) 
			if (is.end){
				n.path <- n.path + 1
				path.list[[n.path]] <- rev(unlist(as.list(cur.path)))
				print(path.list[[n.path]])
				readline()
			}
			if (is.end | is.visited){
				tmp <- pop(cur.path)
				cur.node <- node.prev[cur.node]
				next
			}
			
			# move to the current daughter node and record the current node in node.prev
			cur.daughter <- adj.list[[cur.node]][node.next[cur.node]]
			node.prev[cur.daughter] <- cur.node
			adj.list.visited[[cur.node]][node.next[cur.node]] <- TRUE
			node.next[cur.node] <- node.next[cur.node] + 1
			cur.node <- cur.daughter
			push(cur.path, cur.node)
		}
	}
}

