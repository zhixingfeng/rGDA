encode.snp <- function(x)
{
	base.shift <- integer(length(x$read_base))
	base.shift[x$read_base=='A'] <- 0
	base.shift[x$read_base=='C'] <- 1
	base.shift[x$read_base=='G'] <- 2
	base.shift[x$read_base=='T'] <- 3
	
	x.code <- 4*(x$ref_pos - 1) + base.shift			
	x.code
}

