library(xgboost)
library(Matrix)
context.encode <- function(context)
{
	context.c <- strsplit(context, "")
	
	#context.code <- matrix(0, nrow = length(context.c), ncol = 4*length(rle(context.c[[1]])$values))
	context.code <- Matrix(0, nrow = length(context.c), ncol = 4*length(rle(context.c[[1]])$values))

	for(i in 1:length(context.c)){
		if ((i+1)%%10000 == 0) print (i+1)
		cur.context.count <- rle(context.c[[i]])
		len <- length(cur.context.count$values)
		
		shift <- integer(len)
		shift[cur.context.count$values == 'A'] <- 0
		shift[cur.context.count$values == 'C'] <- 1
		shift[cur.context.count$values == 'G'] <- 2
		shift[cur.context.count$values == 'T'] <- 3
		
		context.code[i, 4*(0:(len-1)) + shift + 1] <- cur.context.count$lengths
	}
	context.code
}

context.effect.encode <- function(context.effect.data, is.train=FALSE)
{
	# remove context containing uncertain bases
	context.raw <- apply(context.effect.data, 1, function(x) paste(x[1], x[2], sep=""))
	idx.rm <- which(grepl("[R|Y|S|W|K|M|B|D|H|V|N]", context.raw))
	
	context <- context.raw
	context.effect.data.valid <- context.effect.data
	if (length(idx.rm) > 0){
		context <- context.raw[-idx.rm]
		context.effect.data.valid <- context.effect.data[-idx.rm]
	}

	# encode context
	context.code <- context.encode(context)
	
	# build 4 models for A, C, G, T respectively
	len.left <- sapply(context.effect.data.valid[,1], nchar)
	ref.base <- substr(context.effect.data.valid[,1], len.left, len.left)
	
	# model A
	idx.A <- which(ref.base != "A")	
	freq.A <- context.effect.data.valid[idx.A, 3] / context.effect.data.valid[idx.A, 7]
	context.code.A <- context.code[idx.A, ]
	#model.A.cv <- xgb.cv(data = context.code.A, label = freq.A, nrounds = 100, nfold = 5, max_depth = 6, objective = "reg:linear")	
	
	# model C
	idx.C <- which(ref.base != "C")
        freq.C <- context.effect.data.valid[idx.C, 3] / context.effect.data.valid[idx.C, 7]
        context.code.C <- context.code[idx.C, ]

	# model G
	idx.G <- which(ref.base != "G")
        freq.G <- context.effect.data.valid[idx.G, 3] / context.effect.data.valid[idx.G, 7]
        context.code.G <- context.code[idx.G, ]

	# model T	
	idx.T <- which(ref.base != "T")
        freq.T <- context.effect.data.valid[idx.T, 3] / context.effect.data.valid[idx.T, 7]
        context.code.T <- context.code[idx.T, ]
	
	if (is.train){
	}else{
		list(context.code = context.code, idx.A = idx.A, idx.C = idx.C, idx.G = idx.G, idx.T = idx.T,
			freq.A = freq.A, freq.C = freq.C, freq.G = freq.G, freq.T = freq.T)
	}
}

context.effect.learn.cv <- function(context.effect.encode.obj, nrounds = 400, eta = 0.01, nfold = 5, max_depth = 1, subsample = 0.5, min_child_weight = 5, objective = "reg:linear")
{
	cv.A <- xgboost.cv(dat = context.effect.encode.obj$context.code[context.effect.encode.obj$idx.A,],
                label = context.effect.encode.obj$freq.A, nrounds = nrounds, eta = eta, nfold = nfold, max_depth = max_depth,
                subsample = subsample, min_child_weight = min_child_weight)	
	
	cv.C <- xgboost.cv(dat = context.effect.encode.obj$context.code[context.effect.encode.obj$idx.C,],
                label = context.effect.encode.obj$freq.C, nrounds = nrounds, eta = eta, nfold = nfold, max_depth = max_depth,
                subsample = subsample, min_child_weight = min_child_weight)

	cv.G <- xgboost.cv(dat = context.effect.encode.obj$context.code[context.effect.encode.obj$idx.G,],
                label = context.effect.encode.obj$freq.G, nrounds = nrounds, eta = eta, nfold = nfold, max_depth = max_depth,
                subsample = subsample, min_child_weight = min_child_weight)

	cv.T <- xgboost.cv(dat = context.effect.encode.obj$context.code[context.effect.encode.obj$idx.T,],
                label = context.effect.encode.obj$freq.T, nrounds = nrounds, eta = eta, nfold = nfold, max_depth = max_depth,
                subsample = subsample, min_child_weight = min_child_weight)

	list(cv.A = cv.A, cv.C = cv.C, cv.G = cv.G, cv.T = cv.T)
}

#nrounds = 400, eta = 0.01, nfold = 5, max_depth = 1, subsample = 0.5, min_child_weight = 5, objective = "reg:linear"







