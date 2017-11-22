correct.reads <- function(encode.data, m5.data, var.loci, out.file)
{
	encode.data.correct <- lapply(encode.data, function(x, t) intersect(x,t), t = var.loci)	
	encode.data.idx <- lapply(encode.data.correct, function(x,t) match(x,t), t = var.loci)
	encode.data.mat <- matrix(0, nrow=length(encode.data.idx), ncol=length(var.loci))	
	encode.start <- m5.data$tStart*4
	encode.end <- m5.data$tEnd*4 + 3
	for (i in 1:nrow(encode.data.mat)){
		encode.data.mat[i, encode.data.idx[[i]]] <- 1
		encode.data.mat[i, var.loci < encode.start[i] | var.loci > encode.end[i]] <- -1
		
	}
	png(out.file, height=nrow(encode.data.mat)/10, width = ncol(encode.data.mat))		
		image(t(encode.data.mat), col=c("black","white", "blue"), bty='n')
	dev.off()
}






