load.snp.code <- function(snp.mat.file)
{
	load(snp.mat.file)
	# encode snp_mat
	snp.code <- list()
	for (i in 1:nrow(snp.mat)){
        	code.A <- 4 * (which(snp.mat[i,]=='A') - 1) + 0
	        code.C <- 4 * (which(snp.mat[i,]=='C') - 1) + 1
        	code.G <- 4 * (which(snp.mat[i,]=='G') - 1) + 2
	        code.T <- 4 * (which(snp.mat[i,]=='T') - 1) + 3

        	snp.code[[i]] <- sort(c(code.A, code.C, code.G, code.T))
	}

	snp.code
}

eval.centroid <- function(centroids, snp.code, is.fdr = FALSE)
{
	acc <- numeric(length(centroids))
	group.id <- numeric(length(centroids))
	for (i in 1:length(centroids)){
		if (i %% 100 == 0)
			print(i)
        	if (length(centroids[[i]])==0){
                	acc[i] <- NaN
                	next
        	}
        	for (j in 1:length(snp.code)){
                	cur.snp.code <- snp.code[[j]][snp.code[[j]] >= min(centroids[[i]]) & snp.code[[j]] <= max(centroids[[i]])]
                	if (is.fdr){
				cur.acc <- length(intersect(centroids[[i]], cur.snp.code)) / length(centroids[[i]])
			}else{
				cur.acc <- length(intersect(centroids[[i]], cur.snp.code)) / length(union(centroids[[i]], cur.snp.code))
			}
			if (cur.acc > acc[i]){
                        	acc[i] <- cur.acc
                        	group.id[i] <- j
                	}
        	}
	}
	list(acc = acc, group.id = group.id)
}

eval.consensus <- function(cons, snp.code.raw, min.read.count = 10, is.fdr = FALSE)
{
	cons.exl <- which(cons$read.count < min.read.count) - 1
	cons.code <- setdiff(cons$cons.seq, cons.exl)
	
	#snp.code <- setdiff(snp.code.raw, cons.exl)
	snp.code <- lapply(snp.code.raw, function(x,y) setdiff(x,y), y= cons.exl)	

	if (length(cons.code)==0){
        	return(list(acc = NaN, group.id = NaN))
        }

	acc <- 0 
	group.id <- NaN
	for (j in 1:length(snp.code)){
		cur.snp.code <- snp.code[[j]][snp.code[[j]] >= min(cons.code) & snp.code[[j]] <= max(cons.code)]
		if (is.fdr){
                	cur.acc <- length(intersect(cons.code, cur.snp.code)) / length(cons.code)
                }else{
                        cur.acc <- length(intersect(cons.code, cur.snp.code)) / length(union(cons.code, cur.snp.code))
                }
		if (cur.acc >= acc){
			acc <- cur.acc
			group.id <- j
		}
	}
	list(acc = acc, group.id = group.id)
}





