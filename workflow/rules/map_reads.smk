# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: 

###########################
# default

rule annotate_genomes:
    input:
        concatenate_fasta
        fa_index
        bwa_map
        sort_index_bam
        merge_separate_bam
        picard_mark_dup
        gatk_realign_indels
    output:
        touch("status/annotate_genomes.done")


###########
# index separate genomes for the last steps (mark duplicates, etc.)
rule fa_index:
    input:
        drep_genomes
    output:
        'concat_genomes.fasta.fai'
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_index.log")
    threads:
    params:
    conda:
    message:
        "Indexing {wildcards.mag}"
    shell:
        "(date && gatk-launch CreateSequenceDictionary -R {mag}.fa && date) &> {log}"

rule concatenate_fasta:
    input:
        drep_genomes
    output:
    conda:
    shell:
        "cat .. > concat_genomes.fasta"

# index the concat fasta file for mapping
rule fa_index:
    input:
        'concat_genomes.fasta'
    output:
        'concat_genomes.fasta.fai'
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
        ref_genome='concat_genomes.fasta',
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
        "Mapping {wildcards.sample} onto the concatenated fasta file"
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

rule merge_separate_bam:
    input:
        expand(os.path.join(RESULTS_DIR, "mapped_reads/concat_genomes_{sample}.sorted.bam"), sample=SAMPLES)
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
        """(date && 
        # EXAMPLE TO SEPARATE CONTIGS 
        samtools merge merged.bam <(samtools view -b SE.bam chr1) & (samtools view -b SE.bam chr2) 
        
        samtools merge -rh {output} {input} && date) &> {log}"""

rule picard_mark_dup:
    input:
        drep_genomes
    output:
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    shell:
        "java -jar picard.jar MarkDuplicates I={input}.bam O={input}_marked_duplicates.bam M={input}_marked_dup_metrics.txt"

rule gatk_realign_indels:
    input:
        drep_genomes
    output:
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    shell:
        "java –jar GenomeAnalysisTK.jar –T RealignerTargetCreator –R {input}.fa –I {input}_marked_duplicates.bam –o {input}_realigner.intervals"
        "java –jar GenomeAnalysisTK.jar –T IndelRealigner –R {input}.fa –I {input}_marked_duplicates.bam –targetIntervals {input}_realigner.intervals –o {input}_marked_realigned.bam"






