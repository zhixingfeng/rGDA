# calculate read pairs sharing lots of variants
find_read_pairs <- function(encode.data, m5.data, min.overlap.var = 1, is.undirected=TRUE)
{
	###--------- check if encode.data and m5.data match --------###
        if (length(encode.data)!=nrow(m5.data)){
                stop('encode.data and m5.data do not match.')
        }

        ###--------- remove reads contained in another read ------###
        read.pair.mat <- matrix(0, length(encode.data), length(encode.data))
        for (i in 1:length(encode.data)){
                print(i)
                for (j in 1:length(encode.data)){
                        if (i==j) next
                        if (m5.data$tStart[i] > m5.data$tEnd[j] | m5.data$tEnd[i] < m5.data$tStart[j])
                                next
                        n.common.var <- length(intersect(encode.data[[i]], encode.data[[j]]))
                        n.var.overlap <- sum(encode.data[[i]]>=4*m5.data$tStart[j] & encode.data[[i]]<=4*m5.data$tEnd[j]+3)

                        if (n.var.overlap > min.overlap.var & n.common.var >= ceiling(n.var.overlap/2)){
                                read.pair.mat[i,j] <- 1
				if (is.undirected){
					read.pair.mat[j,i] <- 1
				}
			}
                }
        }
	read.pair.mat
}


# overlap layout consensus 
olc <- function(encode.data, m5.data, min.overlap.var = 1, min.overlap = 200)
{
	###--------- check if encode.data and m5.data match --------###
        if (length(encode.data)!=nrow(m5.data)){
                stop('encode.data and m5.data do not match.')
        }

	###--------- remove reads contained in another read ------###
	contain.mat <- matrix(0, length(encode.data), length(encode.data))
	for (i in 1:length(encode.data)){
		print(i)
		for (j in 1:length(encode.data)){
			if (i==j) next
			if (m5.data$tStart[i] < m5.data$tStart[j] | m5.data$tEnd[i] > m5.data$tEnd[j])
				next
			n.common.var <- length(intersect(encode.data[[i]], encode.data[[j]]))
			n.var.overlap <- sum(encode.data[[i]]>=4*m5.data$tStart[j] & encode.data[[i]]<=4*m5.data$tEnd[j]+3)
			
			if (n.var.overlap > 0 & n.common.var >= ceiling(n.var.overlap/2))			
				contain.mat[i,j] <- 1
		}
	}
	is.contained <- apply(contain.mat, 1, function(x) any(x==1))
	idx <- which(!is.contained)
	#save(contain.mat, file = './contain.mat.Rdata')	

	###--------- construct overlap graph ------###	
	overlap.mat <- matrix(0, length(idx), length(idx))
	for (i in 1:length(idx)){
		for (j in 1:length(idx)){
			if (i==j) next
			if (m5.data$tEnd[idx[i]] > m5.data$tEnd[idx[j]])
				next
			n.common.var <- length(intersect(encode.data[[idx[i]]], encode.data[[idx[j]]]))
                        n.var.overlap <- sum(encode.data[[idx[i]]]>=4*m5.data$tStart[idx[j]] &
						 encode.data[[idx[i]]]<=4*m5.data$tEnd[idx[j]]+3)
			if (n.common.var >= ceiling(n.var.overlap/2) &
				n.var.overlap >= min.overlap.var & 
				m5.data$tEnd[idx[i]] - m5.data$tStart[idx[j]] >= min.overlap)
				overlap.mat[i,j] <- 1
			
		}
	}
	list(contain.mat=contain.mat, idx=idx, overlap.mat=overlap.mat)	
}

