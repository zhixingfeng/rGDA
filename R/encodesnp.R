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

decode.snv <- function(x)
{
	locus <- floor(x/4)
	shift <- x%%4
	base <- character(length(x))	
	base[shift==0] <- 'A'
	base[shift==1] <- 'C'
	base[shift==2] <- 'G'
	base[shift==3] <- 'T'
	list(base=base, locus=locus)
}

