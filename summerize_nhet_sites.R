infile <- snakemake@input [["infile"]]
pergene_file <- snakemake@output [["pergene"]]
persample_file<- snakemake@output [["persample"]]


inf<- read.table (infile, head=T)

pergene <- apply (inf [, -c(1)],1,mean)
persample <- apply (inf [, -c(1)],2,mean)

pergene <- data.frame (gene=inf[,1], mean=pergene)
persample <- data.frame (sample = colnames (inf) [-1], mean=persample)

write.table (pergene, pergene_file, col.names=T,row.names=F,quote=F,sep="\t")
write.table (persample, persample_file, col.names=T,row.names=F,quote=F,sep="\t")
