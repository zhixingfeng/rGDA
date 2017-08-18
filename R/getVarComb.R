getVarComb <- function(dforest.data, minfreq=0.65)
{
	dforest.data.sig <- dforest.data[dforest.data[,3]>=minfreq,]
}

