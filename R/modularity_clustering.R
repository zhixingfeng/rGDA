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
mod.cluster <- function(A)
{
	k <- apply(A,1,sum)
        m <- sum(k) / 2
        A.exp <- outer(k,k) / (2*m)
        B <- A - A.exp
	rl.eigs <- eigs_sym(B, 1, "LA")
        membership <- list(gp1=which(rl.eigs$vectors>=0), gp2=which(rl.eigs$vectors<0))

} 

mod.cluster.sub <- function(B.g)
{
	B.delta <- matrix(0, nrow(B.g), ncol(B.g))
	diag(B.delta) <- apply(B.g, 1, sum)
	B.g.delta <- B.g - B.delta
	rl.eigs <- eigs_sym(B.g.delta, 1, "LA")	
	membership <- list(gp1=which(rl.eigs$vectors>=0), gp2=which(rl.eigs$vectors<0))
	list(rl.eigs=rl.eigs, membership=membership, B=B.g)	
}


mod.cluster.core <- function(A)
{
	k <- apply(A,1,sum)
	m <- sum(k) / 2
	A.exp <- outer(k,k) / (2*m)
	#A.exp[A==0] <- 0
	B <- A - A.exp
	rl.eigs <- eigs_sym(B, 1, "LA")
	membership <- list(gp1=which(rl.eigs$vectors>=0), gp2=which(rl.eigs$vectors<0))	
	list(rl.eigs=rl.eigs, membership=membership, B=B)
}

