xgboost.cv <- function(dat, label, nrounds = 400, eta = 0.01, nfold = 5, max_depth = 1, subsample = 0.5, min_child_weight = 5, objective = "reg:linear")
{
	if (nrow(dat) != length(label)) stop ("nrow(dat) != length(label)")
	group.id <- sample(1:nfold, nrow(dat), replace = TRUE)
	pred <- numeric(nrow(dat))
	models <- list()
	for (i in 1:nfold){
		idx <- which(group.id == i)
		dat.test <- xgb.DMatrix(data = dat[idx,], label = label[idx])
		dat.train <- xgb.DMatrix(data = dat[-idx,], label = label[-idx])
		models[[i]] <- xgboost(data = dat.train, nrounds = nrounds, eta = eta, max_depth = max_depth, subsample = subsample, min_child_weight = min_child_weight, objective = objective)
		pred[idx] <- predict(models[[i]], dat.test)
	}
	list(pred = pred, models = models)
}

