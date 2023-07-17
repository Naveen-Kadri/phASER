import os
from collections import defaultdict
infiles = snakemake.input
outfile = snakemake.output.outfile
out = open(outfile, "w")

nvariants = defaultdict(dict)
samples = dict()
nfiles = len(infiles)
for fn, infile in enumerate(infiles):
    print(f'{fn} / {nfiles} - > {fn/nfiles:%}')
    sample = os.path.basename(infile)[:19]
    samples[sample] = samples.get(sample, 0) + 1
    with open(infile, "rt") as inf:
        for lnum, line in enumerate(inf):
            spl = line.rstrip().split()
            if lnum > 0:
                gene, count = spl[3], int(spl[8])
                nvariants[gene][sample] = count

header = "\t".join(["gene"] + list(samples.keys()))
out.write(f'{header}\n')

for gene in nvariants:
    tw = []
    for sample in samples:
        tw.append(str(nvariants[gene].get(sample, 0)))
    tw = "\t".join(tw)
    out.write(f'{gene}\t{tw}\n')
out.close()
