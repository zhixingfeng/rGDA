conv_encode_to_mat <- function(encode.data)
{
	var.set <- sort(unique(unlist(encode.data)))
	encode.mat <- matrix(0, nrow = length(encode.data), ncol = length(var.set))
	if (nrow(encode.mat)==0 | ncol(encode.mat)==0)
		return(NULL)
	encode.data.match <- lapply(encode.data, function(x,t) match(x,t), t = var.set)
	for (i in 1:length(encode.data.match)){
		encode.mat[i,encode.data.match[[i]]] <- 1
	}
	encode.mat
}


