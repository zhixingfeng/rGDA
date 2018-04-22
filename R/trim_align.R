trim.m5 <- function(encode.data.raw, m5.data.raw)
{
	if (length(encode.data.raw) != nrow(m5.data.raw))
		stop('length(encode.data.raw) != nrow(m5.data.raw)')

	tStart <- floor(sapply(encode.data.raw, function(x) if(length(x)>0) min(x) else NaN ) / 4)
	tEnd <-floor(sapply(encode.data.raw, function(x) if(length(x)>0) max(x) else NaN) / 4)
	idx <- which(!is.na(tStart))
	if (any(tStart[idx] > tEnd[idx]))
		stop('any(tStart > tEnd)')

	m5.data.trim <- m5.data.raw
	for (i in idx){
		m5.data.trim$tStart[i] <- tStart[i]
		m5.data.trim$tEnd[i] <- tEnd[i]
		
		match.pos <- cumsum(s2c(m5.data.raw$tAlignedSeq[i]) != '-')
		trim.start <- m5.data.trim$tStart[i] - m5.data.raw$tStart[i] + 1
		trim.end <- m5.data.trim$tEnd[i] - m5.data.raw$tStart[i] + 1
		
		m5.data.trim$qAlignedSeq[i] <- c2s( s2c(m5.data.raw$qAlignedSeq[i])[match.pos >= trim.start & match.pos <= trim.end] )
		m5.data.trim$tAlignedSeq[i] <- c2s( s2c(m5.data.raw$tAlignedSeq[i])[match.pos >= trim.start & match.pos <= trim.end] )
		m5.data.trim$matchPattern[i] <- c2s( s2c(m5.data.raw$matchPattern[i])[match.pos >= trim.start & match.pos <= trim.end] )
	}
	m5.data.trim
}


