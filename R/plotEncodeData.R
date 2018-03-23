plot.encode.data <- function(encode.data, m5.data, out.file)
{
	var.loci <- sort(unique(unlist(encode.data)))
	dat.mat <- matrix(0, nrow = length(encode.data), ncol = length(var.loci))
	encode.data.idx <- lapply(encode.data, function(x, t) match(x,t), t = var.loci)
	for (i in 1:length(encode.data.idx)){
        	cur.start <- 4*m5.data$tStart[i]
        	cur.end <- 4*m5.data$tEnd[i]
        	dat.mat[i, var.loci>=cur.start & var.loci<=cur.end] <- 1
        	dat.mat[i, encode.data.idx[[i]]] <- 2
	}
	png(out.file)
		image(t(dat.mat), col=c('black','white','red'))
	dev.off()

}


