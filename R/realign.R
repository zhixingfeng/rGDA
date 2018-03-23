const.haplotype <- function(centroids, ref.seq)
{
        haplo.seq <- list()
	for (i in 1:length(centroids)){
                centroids.decode <- decode.snv(centroids[[i]])
                haplo.seq[[i]] <- ref.seq[[1]]
                haplo.seq[[i]][centroids.decode$locus + 1] <- centroids.decode$base
        }
	names(haplo.seq) <- paste('centroid_', 1:length(centroids), sep='')
	haplo.seq
}

#const.haplotype <- function(centroids, ref.seq, out.dir)
#{
#	system(paste('mkdir -p', out.dir))
#	for (i in 1:length(centroids)){
#		centroids.decode <- decode.snv(centroids[[i]])
#		haplo.seq <- ref.seq
#		haplo.seq[[1]][centroids.decode$locus + 1] <- centroids.decode$base 
#		haplo.file <- paste(out.dir, '/centroid_', i, '.fa', sep='')
#		write.fasta(haplo.seq, paste('centroid_', i, sep=''), file.out = haplo.file)
#	}
#}



