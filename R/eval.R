eval.detection <- function(detected.code, true.code, exclude.loci = integer())
{
	exclude.code <- c(4*exclude.loci, 4*exclude.loci+1, 4*exclude.loci+2, 4*exclude.loci+3)

	detected.code.ft <- setdiff(detected.code, exclude.code)
	true.code.ft <- setdiff(true.code, exclude.code)

	n.common <- length(intersect(detected.code.ft, true.code.ft))
	
}

split.ann.by.group <- function(rl.eval, true.encode)
{
	lapply(1:length(true.encode), function(x, y)  which(sapply(y, function(t,k) any(t==k), k = x)), y = rl.eval$group.id)
}


eval.ann.heatmap <- function(ann.data, true.encode, var.data, pdf.file = "")
{
	rl.eval <- eval.ann(ann.data, true.encode)
	ann.id.gp <- split.ann.by.group(rl.eval, true.encode)

	all.loci <- sort(unique(c(floor(unlist(true.encode)/4), var.data$locus)))

	# merge all contigs corresponding to each true contig
	ann.merged <- list()
	for (i in 1:length(ann.id.gp)){
		ann.merged[[i]] <- list()

		ann.merged[[i]]$tested.loci <- integer(0)
		ann.merged[[i]]$cons.seq <- integer(0)
		ann.merged[[i]]$start <- -1
		ann.merged[[i]]$end <- -1
		
		if (length(ann.id.gp[[i]]) == 0)
			next

		ann.merged[[i]]$tested.loci <- sort(unique(unlist(ann.data$tested_loci[ann.id.gp[[i]]])))
		ann.merged[[i]]$cons.seq <- sort(unique(unlist(ann.data$cons_seq[ann.id.gp[[i]]])))	
	}
	
	if (length(ann.merged) != length(true.encode))
		stop('length(ann.merged) != length(true.encode)')
	
	if (pdf.file == ""){
		return(list(rl.eval = rl.eval, ann.id.gp = ann.id.gp, all.loci = all.loci, ann.merged = ann.merged))
	}

	# draw heatmap
	heatmap.data <- matrix(-2, nrow = 3*length(true.encode), ncol = length(all.loci))
	rownames(heatmap.data) <- rep("", 3*length(true.encode))
	for (i in 1:length(ann.merged)){
		cur.locus.true.loci <- match(floor(true.encode[[i]] / 4), all.loci)
		cur.locus.contig.loci <- match(floor(ann.merged[[i]]$cons.seq / 4), all.loci)
		cur.locus.contig.tested_loci <- match(ann.merged[[i]]$tested.loci, all.loci)
	
		cur.locus.true.base <- true.encode[[i]] %% 4
		cur.locus.contig.base <- ann.merged[[i]]$cons.seq %% 4

		if (any(is.na(c(cur.locus.true.loci, cur.locus.contig.loci))))		
			stop('any(is.na(c(cur.locus.true, cur.locus.contig)))')

		t <- 2*(i-1)
		heatmap.data[t+1, ] <- 0
		heatmap.data[t+1, cur.locus.true.loci] <- cur.locus.true.base + 1
		
		heatmap.data[t+2, ] <- -1
		heatmap.data[t+2, cur.locus.contig.tested_loci] <- 0
		heatmap.data[t+2, cur.locus.contig.loci] <- cur.locus.contig.base + 1
		
		rownames(heatmap.data)[t+1] <- paste('True contig ',i, sep='')
		rownames(heatmap.data)[t+2] <- paste('Reconstructed contig ',i, sep='')
	}

	pdf(pdf.file, height = 0.15*nrow(heatmap.data))
		#heatmap.2(heatmap.data[1:10,1:100], Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = 'none', col = c('black', 'white', 'green', 'red', 'orange', 'blue'), trace = 'none')
		#image(t(heatmap.data), col = c('black', 'white', 'green', 'red', 'orange', 'blue'))
		#heatmap(heatmap.data, Rowv = NA, Colv = NA, scale="none", col = c('black', 'white', 'green', 'red', 'orange', 'blue'))
		pheatmap(heatmap.data, cluster_rows = FALSE, cluster_cols=FALSE, color = c('white','black', 'grey', 'green', 'red', 'orange', 'blue'))
	dev.off()
	list(rl.eval = rl.eval, ann.id.gp = ann.id.gp, all.loci = all.loci, ann.merged = ann.merged)
}


eval.ann.hamming <- function(ann.data, true.encode, is.var = FALSE)
{
	g.start <- 0
	g.end <- floor(max(unlist(true.encode))/4 + 1)
	m5.2 <- list(tStart = g.start, tEnd = g.end)
	
	var.data <- list(locus = sort(unique( floor(unlist(true.encode)/4) )))
        acc <- rep(0,nrow(ann.data))
        group.id <- rep(NaN,nrow(ann.data))
        for (i in 1:nrow(ann.data)) {
                if (i %% 100 == 0)
                        print(i)
                for (j in 1:length(true.encode)){
                	m5.1 <- list(tStart = ann.data$start[i], tEnd = ann.data$end[i])
			cur.acc <- 1 - dist_hamming(ann.data$cons_seq[[i]], m5.1, true.encode[[j]], m5.2, var.data, is.var)$dist
			if (cur.acc > acc[i]){
                                group.id[i] <- j
                                acc[i] <- cur.acc
                        }	
		}
        }
        list(acc = acc, group.id = group.id)
}

eval.ann <- function(ann.data, true.encode, is.fdr = FALSE, is.trim = TRUE)
{
	acc <- rep(-1,nrow(ann.data))
	group.id <- lapply(1:nrow(ann.data), function(x) integer(0))
	true.encode.trim <- list()
	fp <- list()
	fn <- list()
	for (i in 1:nrow(ann.data)) {
		if (i %% 100 == 0)
			print(i)
		for (j in 1:length(true.encode)){
			if (is.trim){
				cur.true.encode <- intersect(true.encode[[j]], 
				c(4*ann.data$tested_loci[[i]], 4*ann.data$tested_loci[[i]]+1, 4*ann.data$tested_loci[[i]]+2, 4*ann.data$tested_loci[[i]]+3))
			}else{
				cur.true.encode <- true.encode[[j]][floor(true.encode[[j]]/4)>=ann.data$start[i] & floor(true.encode[[j]]/4) <= ann.data$end[i]]
			}
			if (is.fdr){
				cur.acc <- length(intersect(ann.data$cons_seq[[i]], cur.true.encode)) / length(ann.data$cons_seq[[i]])
			}else{
				cur.acc <- length(intersect(ann.data$cons_seq[[i]], cur.true.encode)) / length(union(ann.data$cons_seq[[i]], cur.true.encode))
			}
			if (cur.acc >= acc[i]){
				#group.id[i] <- j
				if (cur.acc > acc[i]){
					group.id[[i]] <- j
					true.encode.trim[[i]] <- cur.true.encode
					fp[[i]] <- setdiff(ann.data$cons_seq[[i]], cur.true.encode)
					fn[[i]] <- setdiff(cur.true.encode, ann.data$cons_seq[[i]])
				}else{
					group.id[[i]] <- c(group.id[[i]], j) 
				}
				acc[i] <- cur.acc
			}
		}
	}
	list(acc = acc, group.id = group.id, true.encode.trim = true.encode.trim, fp = fp, fn = fn)
}

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





