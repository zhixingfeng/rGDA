#library(RSpectra)
#sc <- function(W, method = "raw", k = 2)
#{
#	D <- matrix(0, nrow=nrow(W), ncol=ncol(W))
#	diag(D) <- apply(W, 1, sum)
#	L <- D - W
#	D_inv <- D
#	diag(D_inv) <- 1/diag(D)
#	D_inv_1_2 <- D
#	diag(D_inv_1_2) <- 1/sqrt(diag(D))
#	if (method == "sym")
#		L <- D_inv_1_2 %*% L %*% D_inv_1_2
#	if (method == "rw")
#                L <- D_inv %*% L
#
#	rl <- eigs_sym(L, k, which = "SM")
#	rl 
#	#sort(rl$values)
#}


