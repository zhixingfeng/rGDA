parse.nucmer <- function (nucmer.coords.file) 
{
    nucmer.coords <- read.table(nucmer.coords.file, header = FALSE, 
        as.is = TRUE, sep = "\t")
    names(nucmer.coords) <- c("ref_start", "ref_end", "read_start", 
        "read_end", "ref_map_len", "read_map_len", "Idt", "ref_len", 
        "read_len", "ref_name", "read_name")
    nucmer.coords
}

parse.snp <- function (nucmer.snp.file) 
{
    nucmer.snp <- read.table(nucmer.snp.file, header = FALSE, 
        as.is = TRUE, sep = "\t")
    names(nucmer.snp) <- c("ref_pos", "ref_base", "read_base", 
        "read_pos", "buff", "dist", "na_1", "na_2", "ref_len", 
        "read_len", "na_3", "na_4", "ref_name", "read_name")
    nucmer.snp
}





