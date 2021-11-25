# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment
import os
import glob
import pandas as pd
import subprocess

localrules: 

###########################
rule phylo_analysis:
    input:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/gtotree_out/'), GENUS=GENUS_LIST)
    output:
        touch("status/phylo_analysis.done")

rule run_fastani:
    input:
        GENUS=GENUS_LIST
    output:
        directory(expand(os.path.join(RESULTS_DIR, '{GENUS}/gtotree_out'), GENUS=GENUS_LIST))
    conda:
        os.path.join(ENV_DIR, "fastani.yaml")
    run:
        for i in range(0,len(input.GENUS)):
            paths_file = os.path.join(res_dir, input.GENUS[i], 'paths_list.txt')
            args = ['./FastANI/fastANI','--rl',paths_file,'--ql',paths_file,'-t',32,'-o','fastani.out']

rule run_gtotree:
    input:
        GENUS=GENUS_LIST
    params:
        HMM=HMM_LIST,
        OUT=OUT_LIST,
        RES_DIR=RESULTS_DIR
    output:
        directory(expand(os.path.join(RESULTS_DIR, '{GENUS}/gtotree_out'), GENUS=GENUS_LIST))
    conda:
        os.path.join(ENV_DIR, "gtotree.yaml")
    script:
        os.path.join(SRC_DIR, "run_gtotree.py")
