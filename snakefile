import os
'''
For haplotype based allelel specific expression
phaser requires pysam Cython intervaltree(for gene_ae)
'''
tissues = ['testis', 'epi_h', 'vas_d']
tissues = ['testis', 'epi_h']
chromosomes = range(1, 30)


gtf_file = "/cluster/work/pausch/naveen/ASE/INPUTS/Bos_taurus.ARS-UCD1.2.104.gtf.gz"
vcf = '/cluster/work/pausch/xena/imputation/filtered/{chr}_beagle.vcf.gz'
bam_file = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/{tissue}/RNA_alignments/wasp/wasp_filter/{sample}/{sample}_filter.bam'
OUT_DIR = '/cluster/work/pausch/naveen/ASE/HAP_MAF'
cis_eqtl = '/cluster/work/pausch/xena/eQTL/cis_all/{tissue}/maf01/conditional/FDR_05/conditional_top_variants.txt'
samplefile = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/{tissue}/TPM_matrices/{tissue}_samples_good.txt'


# GET SAMPLE NAMES PER TISSUE
samples = {}
# the second sample added to the latest run! why?
tormv = ['BSWCHEM120151536851', 'BSWCHEM120153327877', 'BSWCHEM120148644897V']
tormv = ['BSWCHEM120148644897V']

for tissue in tissues:
    mysamples = []
    mysamplefile = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/eQTL/{tissue}/TPM_matrices/{tissue}_samples_good.txt'.format(
        tissue=tissue)
    with open(mysamplefile) as inf:
        for line in inf:
            xsample=line.rstrip().split() [0]
            if xsample not in tormv:
                mysamples.append(xsample)
    samples[tissue] = mysamples
for tissue in tissues:
    print(f'Number of samples for {tissue} : {len (samples [tissue])  }')

# TOOLS
# run with python 2.7
PHASER = '/cluster/home/nkadri/ASE/phASER/phaser'
BEAGLE = '/cluster/home/nkadri/PROGRAMS/BEAGLE5.4/beagle.22Jul22.46e.jar '
plink2 = '/cluster/home/nkadri/PROGRAMS/PLINK2/plink2 '
# exts for output files from phaser
output_exts = ['variant_connections.txt', 'allelic_counts.txt', 'haplotypes.txt',
               'haplotypic_counts.txt', 'allele_config.txt', 'vcf.gz', 'vcf.gz.tbi']

rule all:
    input:
        expand(
            OUT_DIR + "/phASER/{tissue}/CHR{chr}/cis_var_effect.txt", chr=chromosomes, tissue=tissues),
        expand(
            OUT_DIR + "/phASER/{tissue}/genome_cis_var_effect_with_eqtl.txt", tissue=tissues),
        expand(OUT_DIR + "/phASER/{tissue}/n_het_sites.txt", tissue=tissues),
        expand(OUT_DIR + "/phASER/{tissue}/n_het_pergene.txt", tissue=tissues),
        expand(
            OUT_DIR + "/phASER/{tissue}/n_het_persample.txt", tissue=tissues)

rule phase_vcf:
    input:
        vcf = vcf
    output:
        vcf = OUT_DIR + '/VCF/chr_{chr}_beagle5.4.vcf.gz'
    params:
        lambda wc, output: output.vcf[:-7],
        lambda wc, resources, threads: int (threads * (resources.mem_mb/1000)),
        window = 60,
        overlap = 10,
        tempdir = '/cluster/work/pausch/temp_scratch/naveen/',
    resources:
        mem_mb = 8000,
        walltime = '04:00'
    threads:
        4
    shell:
        '''
        module load jdk \n
        java -Djava.io.tmpdir={params.tempdir} -Xss5m -Xmx{params[1]}g -jar {BEAGLE} \
        gt={input.vcf} \
        out={params[0]} \
        window={params.window} \
        overlap={params.overlap} \
        nthreads={threads}
        '''

rule phaser:
    input:
        bam_file = bam_file,
        vcf_file = rules.phase_vcf.output.vcf
    output:
        outfiles = expand(
            OUT_DIR + "/phASER/{{tissue}}/CHR{{chr}}/{{sample}}.{ext}", ext=output_exts)
    params:
        lambda wc: wc.sample[:19],
        lambda wc, output: output.outfiles[0].split('.')[0]
    threads:
        6
    resources:
        mem_mb = 2000,
        walltime = '01:00'
    shell:
        '''
        module load gcc/4.8.5 python/2.7.14 samtools/1.12 htslib/1.12 bedtools2/2.30.0\n
        python {PHASER}/phaser/phaser.py \
        --pass_only 1 \
        --vcf {input.vcf_file} \
        --bam {input.bam_file} \
        --paired_end 1 \
        --mapq 60 \
        --baseq 10 \
        --sample {params[0]} \
        --threads {threads} \
        --chr {wildcards.chr} \
        --o {params[1]}
        '''

rule split_gtf:
    input:
        gtf_file = gtf_file
    output:
        out_file = OUT_DIR + "/GENES/chr{chr}.bed"
    script:
        'split_gtf.py'

rule gene_ae:
    '''
    Annotate the ae results -- here gene info is added
    '''
    input:
        hapcounts = OUT_DIR + \
            "/phASER/{tissue}/CHR{chr}/{sample}.haplotypic_counts.txt",
        features = rules.split_gtf.output.out_file
    output:
        outfile = OUT_DIR + "/phASER/{tissue}/CHR{chr}/{sample}_gene_ae.txt"
    resources:
        walltime = '00:10',
        mem_mb = 1000
    shell:
        '''
        module load gcc/4.8.5 python/2.7.14 samtools/1.12 htslib/1.12 bedtools2/2.30.0\n
        python {PHASER}/phaser_gene_ae/phaser_gene_ae.py \
        --haplotypic_counts {input.hapcounts} \
        --features {input.features} \
        --min_haplo_maf 0.05 \
        --o {output.outfile}
        '''


rule count_nhet_sites:
    ''' Count the number of heterozygous sites with reads per sample per gene'''
    input:
        lambda wc: expand(
            OUT_DIR + "/phASER/{{tissue}}/CHR{chr}/{sample}_gene_ae.txt", chr=chromosomes, sample=samples[wc.tissue])
    output:
        outfile = OUT_DIR + "/phASER/{tissue}/n_het_sites.txt"
    resources:
        mem_mb = 6000,
        walltime = '00:20'
    script:
        'count_nhet_sites.py'


rule summerize_nhet_sites:
    input:
        infile = rules.count_nhet_sites.output.outfile
    output:
        pergene = OUT_DIR + "/phASER/{tissue}/n_het_pergene.txt",
        persample = OUT_DIR + "/phASER/{tissue}/n_het_persample.txt"
    script:
        'summerize_nhet_sites.R'


rule all_samples:
    ''' Combine results from all samples into a file --done chromosomewise '''
    input:
        lambda wc: expand(
            OUT_DIR + "/phASER/{{tissue}}/CHR{{chr}}/{sample}_gene_ae.txt", sample=samples[wc.tissue]),
        features = rules.split_gtf.output.out_file
    output:
        outfiles = expand(
            OUT_DIR + "/phASER/{{tissue}}/CHR{{chr}}/all_samples.{ext}", ext=["bed.gz", "gw_phased.bed.gz"]),
    params:
        lambda wc, input: os.path.dirname(input[0]),
        lambda wc, output: output.outfiles[0].split('.')[0]
    resources:
        mem_mb = 4000,
        walltime = '00:30',
    threads:
        1
    shell:
        '''
        module load gcc/4.8.5 python/2.7.14 samtools/1.12 htslib/1.12 bedtools2/2.30.0\n
        python {PHASER}/phaser_pop/phaser_expr_matrix.py \
        --gene_ae_dir {params[0]} \
        --features {input.features} \
        --t {threads} \
        --o {params[1]}
        '''

rule merge_vcfs:
    '''merge the phased vcfs [read based phasing] '''
    input:
        lambda wc: expand(
            OUT_DIR + "/phASER/{{tissue}}/CHR{{chr}}/{sample}.vcf.gz", sample=samples[wc.tissue])
    output:
        outfile = OUT_DIR + "/phASER/{tissue}/CHR{chr}/all_samples.vcf.gz"
    resources:
        mem_mb = 4000,
        walltime = "04:00"
    shell:
        '''
        module load gcc/8.2.0 bcftools \n
        bcftools merge {input} --output {output.outfile} --output-type z
        '''

rule index_vcfs:
    input:
        infile = rules.merge_vcfs.output.outfile
    output:
        outfile = rules.merge_vcfs.output.outfile + '.tbi'
    shell:
        '''
        module load gcc/8.2.0 htslib \n
        tabix -p vcf {input.infile}
        '''

rule make_map:
    ''' File mapping the sample names in vcf and bam file required for cis-var'''
    input:
        samplefile = samplefile
    output:
        outfile = OUT_DIR + "/phASER/{tissue}/mapfile.txt"
    script:
        "make_map.py"

rule make_pairs:
    ''' File with info on gene and variant pairs whose ASE effect to be measured '''
    input:
        infile = cis_eqtl
    output:
        outfile = OUT_DIR + "/phASER/{tissue}/CHR{chr}/pairs.txt"
    script:
        "make_pairs.py"


rule var_eff:
    ''' Do not include the - -chr: does not work without the indexed bed file and format is not right for indexing! '''
    input:
        vcf = rules.merge_vcfs.output.outfile,
        index = rules.index_vcfs.output.outfile,
        bed = OUT_DIR +
        "/phASER/{tissue}/CHR{chr}/all_samples.gw_phased.bed.gz",
        map_file = OUT_DIR + "/phASER/{tissue}/mapfile.txt",
        pair_file = rules.make_pairs.output.outfile
    output:
        outfile = OUT_DIR + "/phASER/{tissue}/CHR{chr}/cis_var_effect.txt"
    resources:
        mem_mb = 4000,
        walltime = "04:00",
        tempdir = '/cluster/work/pausch/naveen/TEMP/{tissue}/chr{chr}'
    params:
        min_cov = 8
    threads:
        4
    shell:
        '''
        module load gcc/4.8.5 python/2.7.14 samtools/1.12 htslib/1.12 bedtools2/2.30.0\n
        python {PHASER}/phaser_pop/nkk_phaser_cis_var.py \
        --bed {input.bed} \
        --vcf {input.vcf} \
        --pair {input.pair_file} \
        --map {input.map_file} \
        --o {output.outfile} \
        --min_cov {params.min_cov} \
        --t {threads} \
        '''

rule combine_eff:
    input:
        infiles = expand(
            OUT_DIR + "/phASER/{{tissue}}/CHR{chr}/cis_var_effect.txt", chr=chromosomes)
    output:
        outfile = OUT_DIR + "/phASER/{tissue}/genome_cis_var_effect.txt"
    script:
        "combine_eff.R"


rule get_gt_counts:
    input:
        vcf = vcf,
    output:
        outfile = OUT_DIR + '/phASER/genotype_counts/chr_{chr}.hardy'
    params:
        lambda wc, output: output.outfile[:-6]
    resources:
        mem_mb = 8000,
        walltime = '01:00'
    shell:
        '''
        {plink2} --cow \
        --vcf {input.vcf} \
        --hardy \
        --out {params[0]}
        '''

rule add_gt_counts:
    input:
        counts = expand(
            OUT_DIR + '/phASER/genotype_counts/chr_{chr}.hardy', chr=chromosomes),
        eqtl = cis_eqtl
    output:
        outfile = OUT_DIR + '/phASER/{tissue}/cis_eqtl_with_gt_counts.txt'
    script:
        'add_gt_counts.py'


rule merge_eqtl_ase:
    ''' all eqls are kept.. '''
    input:
        ase = rules.combine_eff.output.outfile,
        eqtl = rules.add_gt_counts.output.outfile
    output:
        outfile = OUT_DIR + \
            "/phASER/{tissue}/genome_cis_var_effect_with_eqtl.txt"
    script:
        'merge_eqtl_ase.R'
