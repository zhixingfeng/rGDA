# remove isolated nodes (nodes with no edge attached)
rm.isolated.nodes <- function(A)
{
	out.degree <- apply(A, 1, sum)	
	in.degree <- apply(A, 2, sum)
	idx.rm <- which(out.degree + in.degree==0)
	list(A=A[-idx.rm, -idx.rm], idx.rm=idx.rm)
}
# core part of modularity clustering (Newman PNAS 2006)
# input A is adjacent matrix
mod.cluster.core <- function(A)
{
	k <- apply(A,1,sum)
	m <- sum(k) / 2
	A.exp <- outer(k,k) / (2*m)
	#A.exp[A==0] <- 0
	B <- A - A.exp
	rl.eigs <- eigs_sym(B, 1, "LA")
	membership <- list(gp1=which(rl.eigs$vectors>=0), gp2=which(rl.eigs$vectors<0))	
	list(rl.eigs=rl.eigs, membership=membership)
}

