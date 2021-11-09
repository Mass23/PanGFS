# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: 

###########################
# checkpoint to collect representative sequences
rule bin_link:
    input:
        os.path.join(RESULTS_DIR, "genomes/{i}.contigs.fa"),
        os.path.join(RESULTS_DIR, "rep_mags/{i}.contigs.fa")
    output:
        os.path.join(RESULTS_DIR, "bins/{sample}__{i}.contigs.fa")
    wildcard_constraints:
        genomes="|".join(GENOMES)
        #i="(\w\.)+"
    shell:
        "ln -vs {input} {output}"

checkpoint collect_genomes:
    input:
    output:
    shell:
        "ls DREP_OUT_DIR/dereplicated_genomes/*.fasta > genomes.txt"

checkpoint collect_mags:
    input:
    output:
    shell:
        "ls DREP_OUT_DIR/dereplicated_genomes/*.fasta > rep_mags.txt"

def bins(wildcards):
    checkpoint_output = checkpoints.collect_genomes.get(**wildcards).output[0] + checkpoints.collect_mags.get(**wildcards).output[0]
    return expand(os.path.join(RESULTS_DIR, "genomes/{{sample}}__{i}.contigs.fa"),
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.contigs.fa")).i)

rule bin_folder_sample:
    input:
        bins
    wildcard_constraints:
        sample="|".join(SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "logs/{sample}_bin_collection.done"),
    shell:
        "for fname in {input} ; do echo $(basename -s \".contigs.fa\" \"${{fname}}\") ; done > {output}"

rule phylo_analysis:
    input:
        os.path.join(RESULTS_DIR, "logs/phylo_analysis.done")
    output:
    params:
        hmm_gtotree=config['gtotree']['hmm']
    shell:
        "GToTree -a NAME_FOR_THE_RUN -f full_genomes_list.txt -H {params.hmm_gtotree} -j 4 -o PanPolyGFS_tree"


