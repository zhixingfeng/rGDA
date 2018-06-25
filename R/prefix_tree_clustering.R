prefix_tree_cluster <- function(encode.data, m5.data, rl.cons)
{
	encode.data.sorted <- lapply(encode.data, function(x, t) x[order(t[x+1], decreasing = TRUE)], t = rl.cons$var.count)
		
}


