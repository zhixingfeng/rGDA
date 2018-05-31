cutwindow <- function(encode.data, m5.data, window)
{
	if (length(encode.data) != nrow(m5.data))
		stop('length(encode.data) != nrow(m5.data)')
	idx <- which(m5.data$tStart <= window[1] & m5.data$tEnd >= window[2])

	encode.data.win <- encode.data[idx]
	m5.data.win <- m5.data[idx,]
	
	encode.data.win.cut <- lapply(encode.data.win, function(x, s, e) x[x>=4*s & x<=4*e+3], s = window[1], e = window[2])
	m5.data.win.cut <- m5.data.win
	for (i in 1:nrow(m5.data.win.cut)){
		m5.data.win.cut$tStart[i] <- window[1]
		m5.data.win.cut$tEnd[i] <- window[2]
		
		match.pos <- cumsum(s2c(m5.data.win$tAlignedSeq[i]) != '-')
                trim.start <- m5.data.win.cut$tStart[i] - m5.data.win$tStart[i] + 1
                trim.end <- m5.data.win.cut$tEnd[i] - m5.data.win$tStart[i] + 1

                m5.data.win.cut$qAlignedSeq[i] <- c2s( s2c(m5.data.win$qAlignedSeq[i])[match.pos >= trim.start & match.pos <= trim.end] )
                m5.data.win.cut$tAlignedSeq[i] <- c2s( s2c(m5.data.win$tAlignedSeq[i])[match.pos >= trim.start & match.pos <= trim.end] )
                m5.data.win.cut$matchPattern[i] <- c2s( s2c(m5.data.win$matchPattern[i])[match.pos >= trim.start & match.pos <= trim.end] )
	}
	
	list(encode.data = encode.data.win.cut, m5.data = m5.data.win.cut)
}


