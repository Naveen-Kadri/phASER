from collections import defaultdict
counts = snakemake.input.counts
eqtl = snakemake.input.eqtl
outfile = snakemake.output.outfile
out = open(outfile, "w")

res = defaultdict(list)
n = 0
with open(eqtl) as inf:
    for i, line in enumerate(inf):
        spl = line.rstrip().split()
        if i == 0:
            header = spl
            myindex = header.index('var_id')
        else:
            n += 1
            res[spl[myindex]].append(spl)

print(f'number of cis e genes : {n}, {len(res):,}')

header = '\t'.join(header + ["rr", "ra", "aa", "hwe"])
out.write(f'{header}\n')

for infile in counts:
    with open(infile) as inf:
        for i, line in enumerate(inf):
            if i == 0:
                continue
            spl = line.rstrip().split()
            if spl[1] in res:
                myres = res[spl[1]]
                for el in myres:
                    tw = '\t'.join(el + spl[4:7] + [spl[9]])
                    out.write(f'{tw}\n')

out.close()
