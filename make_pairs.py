'''
output format :
gene_id	var_id	var_contig	var_pos	var_ref	var_alt
ENSG00000227232.5	chr1_666028_G_A_b38	chr1	666028	G	A
ENSG00000238009.6	chr1_63671_G_A_b38	chr1	63671	G	A
ENSG00000233750.3	chr1_666028_G_A_b38	chr1	666028	G	A

'''
# tissue='testis'
# infile = f'/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/maf01/conditional/FDR_05/conditional_top_variants.txt'
# outfile ='todel.txt'
# mychr = 19

infile = snakemake.input.infile
outfile = snakemake.output.outfile
mychr = str (snakemake.wildcards.chr)



header = "\t".join ( ['gene_id', 'var_id', 'var_contig', 'var_pos', 'var_ref', 'var_alt'] )
out = open (outfile, 'w')
out.write (f'{header}\n')

with open (infile, "rt") as inf:
    for lnum,line in enumerate(inf):
        spl=line.rstrip().split()
        if lnum==0:
            header=spl
        else:
            info=dict(zip(header, spl))
            var_id= info ['var_id']
            if info['var_chr'] != mychr:
                continue
            if ";" in var_id:
                var_id = var_id.split(';')
                for el in var_id:
                    chr, pos, ref, alt = el.split ("_")
                    tw = "\t".join ([info['phe_id'],  info['var_id'], chr, pos, ref, alt ])
                    out.write (f'{tw}\n')
            else:
                chr, pos, ref, alt = var_id.split("_")
                tw = "\t".join ([info['phe_id'],  info['var_id'], chr, pos, ref, alt ])
                out.write (f'{tw}\n')
out.close()
        
