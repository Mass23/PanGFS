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
MAGS = [line.strip() for line in open("mags_list.txt").readlines()]
SAMPLES = [line.strip() for line in open("samples_list.txt").readlines()]
GENOMES = [line.strip() for line in open("accessions_list.txt").readlines()]

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
