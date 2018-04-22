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

check.m5 <- function(m5.data)
{
	tAlign.len <- sapply(m5.data$tAlignedSeq, function(x) sum(s2c(x)!='-'), USE.NAMES = FALSE)
	read.len <- m5.data$tEnd - m5.data$tStart + 1
	tAlign.len == read.len
}


