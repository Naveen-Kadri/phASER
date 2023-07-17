'''
MAPPING SAMPLE NAMES IN VCF AND THE BED FILE
PHASSER adds "_filter" to names in the .bed files
Also the samplenames have suffixes (_EH and _V) to indicate the tissue 
'''

infile=snakemake.input.samplefile
outfile=snakemake.output.outfile
out= open (outfile, 'w')
header ="\t".join (['vcf_sample', 'bed_sample'])
out.write (f'{header}\n')
with open (infile, "rt") as inf:
    for line in inf:
        sample=line.rstrip().split()[0]
        out.write (f'{sample[0:19]}\t{sample + "_filter"}\n')
out.close()
