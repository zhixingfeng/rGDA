permute.reads <- function(encode.file, pu.file, outfile, seed = 18473)
{
	encode.data <- load.encodefile(encode.file)
	pu.data <- read.table(pu.file, sep = "\t", as.is = TRUE)
}




