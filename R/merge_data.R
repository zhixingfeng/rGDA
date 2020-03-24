merge_encode <- function(m5_fofn_file, encode_fofn_file, out_encode_file)
{
	m5_fofn <- readLines(m5_fofn_file)
	encode_fofn <- readLines(encode_fofn_file)
	
	if (length(m5_fofn) != length(encode_fofn)){
		stop("length(m5_fofn) != length(encode_fofn)")
	}
	
	map_encode <- list()
	for (i in 1:length(m5_fofn)){
		print(m5_fofn[i])
		print(encode_fofn[i])
	
		m5_data <- load.m5file(m5_fofn[i])
		encode_data <- load.encodefile(encode_fofn[i])
		if (nrow(m5_data) != length(encode_data)){
			stop("nrow(m5_data) != length(encode_data)")
		}
		
		for (j in 1:nrow(m5_data)){
			if (is.null(map_encode[[m5_data$qName[j]]])){
				map_encode[[m5_data$qName[j]]] <- encode_data[[j]]
			}else{
				map_encode[[m5_data$qName[j]]] <- sort(union(map_encode[[m5_data$qName[j]]], encode_data[[j]]))
			}
		}
	}
	save.encodefile(map_encode, out_encode_file)
	writeLines(names(map_encode), paste(out_encode_file, ".readname", sep = ""))
}


