library(mummerSV)
getsnpmat <- function(snp.files, ref.name, genome.size)
{
	# load snp data
	snp.list <- list()
	snp.locus <- list()
	snp.base <- list()
	snp.ref <- list()
	snp.mat <- matrix('N', nrow=length(snp.files), ncol=genome.size)
	for (i in 1:length(snp.files)){
		if (file.size(snp.files[i])==0)
			next

		snp.list[[i]] <- parse.snp(snp.files[i])
		snp.locus[[i]] <- snp.list[[i]]$ref_pos[snp.list[[i]]$ref_name==ref.name]
		snp.base[[i]] <- snp.list[[i]]$read_base[snp.list[[i]]$ref_name==ref.name]
		snp.ref[[i]] <- snp.list[[i]]$ref_base[snp.list[[i]]$ref_name==ref.name]
		snp.mat[i, snp.locus[[i]]] <- snp.base[[i]]
	}

	# compare snp
	snp.dist <- matrix(0, nrow=length(snp.locus), ncol=length(snp.locus))
	for (i in 1:length(snp.locus)){
		for (j in 1:length(snp.locus)){
			snp.dist[i,j] <- length(setdiff(snp.locus[[i]], snp.locus[[j]])) + length(setdiff(snp.locus[[j]], snp.locus[[i]]))
		}
	}

	list(locus=snp.locus, base=snp.base, ref=snp.ref, dist=snp.dist, mat=snp.mat)
}

