getReadsByID <- function(m5.file, test.file, out.file)
{
	m5.data <- readLines(m5.file)
	test.data <- read.table(test.file, header=FALSE, sep = '\t', as.is = TRUE)
	idx <- unique(test.data[,1]) + 1
	writeLines(m5.data[idx], out.file)
}

