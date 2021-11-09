# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: 

###########################
# default        

rule gtdbtk:
    input:
        os.path.join(RESULTS_DIR, "bins/bin_collection.done")
    output:
        directory(os.path.join(RESULTS_DIR, "gtdbtk_output"))
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    conda:
        os.path.join(ENV_DIR, "gtdbtk.yaml")
    params:
        config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    message:
        "Running GTDB toolkit on MAGs"
    shell:
        "(date && export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fa --genome_dir $(dirname {input}) --out_dir {output} && date) &> >(tee {log})"

rule select_mags:
        
rule create_mags_dir:
    output:
        os.path.join(RESULTS_DIR, "mags_list.txt")
    shell:
        "ls mags/* > mags_list.txt"
        
rule drep_mags:
    input:
        "status/dn_genomes.txt"
    output:
        os.path.join(RESULTS_DIR, "")
    shell:
        "cat genomes_list.txt mags_list.txt > merged_genomes_mags_list.txt"
        "dRep dereplicate PanPolyGFS -p CORES -g merged_genomes_mags_list.txt -sa 0.95 -comp 90 -con 10"


