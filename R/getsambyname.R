getsambyname <- function(name.list, bamfile, outsamfile)
{
	tmp.file <- paste(outsamfile,'.tmp', sep = '')	
	write.table(name.list, quote=FALSE, row.names=FALSE, col.names=FALSE, file = tmp.file)
	
	cmd <- paste("samtools view -H", bamfile, "| sort -u >>", outsamfile)
	print(cmd);system(cmd)
	
	cmd <- paste("fgrep -w -f", tmp.file, bamfile, ">>", outsamfile)
	print(cmd);system(cmd)
	
	cmd <- paste("sam2bam", outsamfile, tmp.file)
	print(cmd);system(cmd)
}

