# Pipeline for mapping reads to bins and then identify SNPs
#
# Example call: snakemake -s Snakefile --configfile config.yaml --use-conda --cores 1 -rpn

##############################
# MODULES
import os, re
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO

##############################
# CONFIG
# can be overwritten by using --configfile <path to config> when calling snakemake
configfile:"config.yaml"

##############################
# Relevant directories
DATA_DIR = config["work_dir"]
RESULTS_DIR = config["results_dir"]
ENV_DIR = config["env_dir"]
SRC_DIR = config["scripts_dir"]
BIN_DIR = config["bin_dir"]
FASTQ_DIR = config["fastq_dir"]

##############################
# Input
MAG = [line.strip() for line in open("polaromonas_mag_list").readlines()]
SAMPLES = [line.strip() for line in open("sample_list").readlines()]
SED_GI = [line.strip() for line in open("sed_gi_list").readlines()]

##############################
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]


##############################
rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.bam"), mag=SED_GI, sample=SAMPLES),
#        expand(os.path.join(RESULTS_DIR, "merged_bam/{mag}_merged.bam"), mag=MAG),
#        expand(os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf.gz"), mag=MAG),
#        expand(os.path.join(RESULTS_DIR, "coverage/{mag}_depth.cov"), mag=SED_GI),
        expand(os.path.join(RESULTS_DIR, "mapped_reads/{sed_gi}_{sample}.genome.txt"), sed_gi=SED_GI, sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "sed_gi_coverage/{sed_gi}_{sample}_depth.cov"), sed_gi=SED_GI, sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "sed_gi_vcf/{sed_gi}_{sample}_filtered.bcf.gz"), sed_gi=SED_GI, sample=SAMPLES)


###########
# Mapping #
rule fa_index:
    input:
        ref_genome=os.path.join(BIN_DIR, "{mag}.fa"),
    output:
        expand(os.path.join(BIN_DIR, "{{mag}}.{ext}"), mag=MAG, ext=BWA_IDX_EXT)
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_index.log")
    threads:
        config["bwa"]["threads"]
    params:
        idx_prefix=lambda wildcards, output: os.path.splitext(output[0])[0]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Indexing {wildcards.mag}"
    shell:
        "(date && bwa index {input} -p {params.idx_prefix} && date) &> {log}"

rule bwa_map:
    input:
        ref_genome=os.path.join(BIN_DIR, "{mag}.fa"),
        r1=os.path.join(FASTQ_DIR, "{sample}_R1.fastp.fastq.gz"),
        r2=os.path.join(FASTQ_DIR, "{sample}_R2.fastp.fastq.gz"),
        idx=expand(os.path.join(BIN_DIR, "{mag}.{ext}"), mag=MAG, ext=BWA_IDX_EXT)
    output:
        os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.bam")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_{sample}.log")
    threads:
        config["bwa"]["threads"]
    params:
        idx_prefix=lambda wildcards, input: os.path.splitext(input.idx[0])[0]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    wildcard_constraints:
        mag="|".join(MAG)
    message:
        "Mapping {wildcards.sample} onto {wildcards.mag}"
    shell:
        """(date && bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' $(echo {input.ref_genome} | sed 's/.fa//g') {input.r1} {input.r2} | samtools view -Sb -F 4 - > {output} && date) &> {log}"""

rule sort_index_bam:
    input:
        os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.bam")
    output:
        sorted=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.sorted.bam"),
        index=os.path.join(RESULTS_DIR, "mapped_reads/{mag}_{sample}.sorted.bam.bai")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_{sample}_bam_index.log")
    threads:
        config["bwa"]["sort"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    params:
        chunk_size=config["bwa"]["sort"]["chunk_size"]
    message:
        "Sorting and indexing bam files for {wildcards.mag} and {wildcards.sample}"
    shell:
        "(date && samtools sort --threads {threads} -m {params.chunk_size} {input} > {output.sorted} && "
        "samtools index {output.sorted} && date) &> {log}"

rule merge_bam:
    input:
        expand(os.path.join(RESULTS_DIR, "mapped_reads/{{mag}}_{sample}.sorted.bam"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "merged_bam/{mag}_merged.bam")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_merged_bam.log")
    threads:
        config["bwa"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Merging the bam files for all samples for {wildcards.mag}"
    shell:
        "(date && samtools merge -rh {output} {input} && date) &> {log}"

rule bam2vcf:
    input:
        ref=os.path.join(BIN_DIR, "{mag}.fa"),
        bam=rules.merge_bam.output
    output:
        gz=os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_bam2vcf.log")
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    params:
        calls=temp(os.path.join(RESULTS_DIR, "vcf/{mag}_calls.bcf")),
        filtered=temp(os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf"))
    message:
        "Calling SNPs on {wildcards.mag}"
    shell:
        "(date && mkdir -p $(dirname {output.gz}) && bcftools mpileup -a FORMAT/AD,FORMAT/DP --max-depth 20000 -f {input.ref} {input.bam} | bcftools call -mv -Ob -o {params.calls} && "
        "bcftools view -i '%QUAL>=30' {params.calls} > {params.filtered} && "
        "bgzip {params.filtered} > {output.gz} && date) &> {log}"

rule sed_GI_bam2vcf:
    input:
        ref=os.path.join(BIN_DIR, "{sed_gi}.fa"),
        bam=os.path.join(RESULTS_DIR, "mapped_reads/{sed_gi}_{sample}.sorted.bam")
    output:
        gz=os.path.join(RESULTS_DIR, "sed_gi_vcf/{sed_gi}_{sample}_filtered.bcf.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/{sed_gi}_{sample}_sed_GI_vcf.log")
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    threads:
        config["bwa"]["vcf"]["threads"]
    wildcard_constraints:
        sed_gi="|".join(SED_GI)
    params:
        calls=temp(os.path.join(RESULTS_DIR, "sed_gi_vcf/{sed_gi}_{sample}_calls.bcf")),
        filtered=temp(os.path.join(RESULTS_DIR, "sed_gi_vcf/{sed_gi}_{sample}_filtered.bcf"))
    message:
        "Calling SNPs on {wildcards.sed_gi} and {wildcards.sample}"
    shell:
        "(date && mkdir -p $(dirname {output.gz}) && bcftools mpileup -a FORMAT/AD,FORMAT/DP --threads {threads} --max-depth 20000 -f {input.ref} {input.bam} | bcftools call -mv -Ob -o {params.calls} && "
        "bcftools view -i '%QUAL>=30' {params.calls} > {params.filtered} && "
        "bgzip -f {params.filtered} > {output.gz} && date) &> {log}"


###########
# Depth
rule genome_file:
    input:
        os.path.join(RESULTS_DIR, "mapped_reads/{sed_gi}_{sample}.sorted.bam")
    output:
        os.path.join(RESULTS_DIR, "mapped_reads/{sed_gi}_{sample}.genome.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/{sed_gi}_{sample}_genome.log")
    threads:
        config["bwa"]["sort"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Creating GENOME file for bedtools for {wildcards.sed_gi} and {wildcards.sample}"
    shell:
        "(date && samtools view -H {input} | grep @SQ | sed 's/@SQ\\tSN:\\|LN://g' > {output} && date) &> {log}"

rule depth:
    input:
#        ref=os.path.join(BIN_DIR, "{mag}.fa"),
#        bam=expand(os.path.join(RESULTS_DIR, "mapped_reads/{{mag}}_{sample}.bam"), sample=SAMPLES)
        gen_file=os.path.join(RESULTS_DIR, "mapped_reads/{sed_gi}_{sample}.genome.txt"),
        bam=os.path.join(RESULTS_DIR, "mapped_reads/{sed_gi}_{sample}.sorted.bam")
    output:
#        wind=temp(os.path.join(RESULTS_DIR, "coverage/{mag}_SlidingWindows1kb.bed")),
#        depth=os.path.join(RESULTS_DIR, "coverage/{mag}_depth.cov")
        wind=temp(os.path.join(RESULTS_DIR, "sed_gi_coverage/{sed_gi}_{sample}_SlidingWindows1kb.bed")),
        depth=os.path.join(RESULTS_DIR, "sed_gi_coverage/{sed_gi}_{sample}_depth.cov")
    log:
#        os.path.join(RESULTS_DIR, "logs/{mag}_coverage.log")
        os.path.join(RESULTS_DIR, "logs/{sed_gi}_{sample}_coverage.log")
    threads:
        config["bwa"]["vcf"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bedtools.yaml")
    wildcard_constraints:
        sed_gi="|".join(SED_GI)
    message:
#        "Estimating the coverage for all samples on {wildcards.mag}"
        "Estimating the coverage for {wildcards.sed_gi} and {wildcards.sample}"
    shell:
        "(date && mkdir -p $(dirname {output.depth}) && bedtools makewindows -g {input.gen_file} -w 1000 > {output.wind} && "
        "bedtools multicov -bams {input.bam} -bed {output.wind} > {output.depth} && date) &> {log}"

rule create_dirs:
    input:
    output:
    log:
    threads:
    conda:
    message:
        "Create directories for the next steps."
    shell:


rule annot_genome:
    input:
    output:
    log:
    threads:
    conda:
    message:
        "Annotate the genome with prokka and mantis"
    shell:
        "prokka --prefix {mag} --kingdom Bacteria --genus polaromonas genomes/{mag}.fa"
        "python mantis/ run_mantis -t {mag}/{mag}.faa -o {mag}/ -od polaromonas --mantis_config {input.config} --hmmer_threads {params.cores} --cores {threads} --memory {config[mantis][single_mem]} --kegg_matrix"

rule nuc_div:
    input:
        {genome}
    output:
    run:
        def list_vcfs(genome):
            return(glob.glob('polymorphism/{genome}*.bcf.gz'))

        for genome in genomes:
            samples = list_vcfs(genome)
            print(samples)
        
            for sample1 in samples:
                subprocess.call(' '.join(['bcftools','index',sample]), shell = True)
                args = ['vcftools', '--gzvcf', sample, '--out', sample.replace('.bcf.gz',''), '--window-pi', '1000', '--minQ', '30', '--maf', '0.1']
                subprocess.call(' '.join(args), shell = True)

    log:
    threads:
    conda:
    message:
        "Compute nuc. diversity with VCFtools of the VCF files"
