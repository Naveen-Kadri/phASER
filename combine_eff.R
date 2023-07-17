infiles  <- snakemake@input[["infiles"]]
outfile  <- snakemake@output [["outfile"]]

out  <- data.frame()
for (infile in infiles) {
    inf  <- read.table (infile, head=T,sep="\t")
    out  <- rbind (out, inf)
}

write.table (out, outfile, col.names=T,row.names=F,quote=F,sep="\t")
