validate.data <- function(encode.data, m5.data)
{
	if (length(encode.data) != nrow(m5.data))
		stop('length(encode.data) != nrow(m5.data)')
	is.val <- logical(length(encode.data))
	for (i in 1:length(encode.data)){
		is.val[[i]] <- !any(floor(encode.data[[i]]/4) < m5.data$tStart[i] || floor(encode.data[[i]]/4) > m5.data$tEnd[i])
		if (length(encode.data[[i]]) == 0)
			is.val[[i]] <- TRUE
	}
	is.val
}

