# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: 

###########################
rule phylo_analysis:
    input:
        os.path.join(RESULTS_DIR,"gtotree_output")
    output:
        touch("status/phylo_analysis.done")

rule run_gtotree:
    input:
        GENUS=GENUS_LIST
    output:
        directory(os.path.join(RESULTS_DIR, "gtotree_output"))
    conda:
        os.path.join(ENV_DIR, "gtotree.yaml")
    params:
        hmm_gtotree=config['gtotree']['hmm'],
        threads_gtotree=config['gtotree']['threads']
    shell:
        #"gunzip {RESULTS_DIR}/Genomes/*.fna.gz &&"
        "ls {RESULTS_DIR}/{input.GENUS}/mags/dereplicated_genomes/*.fasta > drep_mags_paths.txt &&"
        "ls {RESULTS_DIR}/{input.GENUS}/Genomes/*.fna > genomes_paths.txt &&"
        "cat drep_mags_paths.txt genomes_paths.txt > genomes_drep_mags_paths.txt &&"
        "GToTree -o {RESULTS_DIR}/{input.GENUS}/gtotree_output -f genomes_drep_mags_paths.txt -H {params.hmm_gtotree} -j {params.threads_gtotree}"
