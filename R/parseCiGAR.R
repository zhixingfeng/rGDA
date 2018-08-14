parse.cigar <- function(cigar)
{
	idx.type <- gregexpr('(M|I|D|N|S|H|P|=|X)',cigar)[[1]]
	if (length(idx.type)>1){
		idx.type.shift <- c(0, idx.type[1:(length(idx.type)-1)])
	}else{
		idx.type.shift <- 0
	}
	
	cigar_num <- as.integer(mapply(function(x,y,t) substr(t,x+1,y-1) , x = idx.type.shift, y = idx.type, t=cigar ))
	cigar_type <- sapply(idx.type, function(x,t) substr(t, x, x) , t = cigar)	
	list(cigar_num = cigar_num, cigar_type = cigar_type)
}

