get.movie <- function(m5.data)
{
	sapply(strsplit(m5.data$qName, '/'), function(x) x[1])
}

