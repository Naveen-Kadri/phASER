ase_file <- snakemake@input[["ase"]]
eqtl_file <- snakemake@input[["eqtl"]]
outfile <- snakemake@output[["outfile"]]



eq  <-read.table (eqtl_file, head=T )
ase <- read.table (ase_file,head=T,sep="\t")

eq$myid  <- paste (eq$phe_id, eq$var_id,sep=":")
ase$myid  <- paste (ase$gene, ase$var_id, sep=":")

##there are duplicates for some reason -- multi threading ?
ase  <- ase [!duplicated (ase$myid),]

##Merge the two if the ase info is not available still keep them or not ? no not necessarya
all  <- merge (eq, ase, by='myid', all.x=T) ; dim (all)

##remove if there are no pvalues... when there is enough info <8 reads for example
## all  <- all [!is.na (all$het_hom_pvalue),]
## all$fdr <- p.adjust (all$het_hom_pvalue, method='fdr')

write.table(all, outfile, col.names=T,row.names=F,quote=F,sep="\t")
