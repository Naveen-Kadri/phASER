import gzip
from collections import defaultdict
from gtf import *

myftype='gene'
gtf_file = snakemake.input.gtf_file
out_file = snakemake.output.out_file
mychr = str(snakemake.wildcards.chr)

out = open (out_file, "w")
infos = ['seqid', 'start', 'end', 'gene_id', 'strand', 'gene_biotype']

##There could be multiple genes with the same start pos.. 
tosort=defaultdict (list)

#class method to set the items and orders of the attributes to be written with the method print_info
gtf.set_infos_to_write(infos)


#header = "\t".join (infos)
#out.write (f'{header}\n')


with gzip.open (gtf_file, 'rt') as inf:
    for line in inf:
        spl =  line.rstrip().split ()
        if spl [0] == mychr and spl [2] == myftype :
            ##For zero-based 
            spl [3] = str (int(spl [3]) - 1 )
            spl [4] = str (int(spl [4]) - 1 )
            mygtf =gtf (line)
            tosort [ int (spl [3]) ].append ( "\t". join (   [ str(el)  for el in   gtf.print_info (mygtf) ]     ))

for mypos in sorted (tosort.keys() ) :
    for el in tosort[mypos]:
        out.write (f'{el}\n')
