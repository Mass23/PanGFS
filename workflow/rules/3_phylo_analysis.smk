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
        expand(os.path.join(RESULTS_DIR, '{GENUS}/fastani_out'), GENUS=GENUS_LIST)
    output:
        touch("status/phylo_analysis.done")

rule run_gtotree:
    input:
        'status/dn_genomes.done',
        'status/gtdbtk.done'
    params:
        GENUS=GENUS_LIST,
        HMM=HMM_LIST,
        OUT=OUT_LIST,
        RES_DIR=RESULTS_DIR
    output:
        directory(expand(os.path.join(RESULTS_DIR, '{GENUS}/gtotree_out'), GENUS=GENUS_LIST))
    conda:
        os.path.join(ENV_DIR, "gtotree.yaml")
    script:
        os.path.join(SRC_DIR, "run_gtotree.py")

rule run_fastani:
    input:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/gtotree_out'), GENUS=GENUS_LIST)
    params:
        GENUS_LIST=GENUS_LIST
    output:
        directory(expand(os.path.join(RESULTS_DIR, '{GENUS}/fastani_out'), GENUS=GENUS_LIST))
    run:
        os.system('source /home/massimo.bourquin/apps/miniconda3/etc/profile.d/conda.sh')
        os.system('conda activate fastani')
        for i in range(0,len(params.GENUS_LIST)):
            paths_file = os.path.join(RESULTS_DIR, input.GENUS_LIST[i], 'paths_list.txt')
            args = ['./FastANI/fastANI','--rl',paths_file,'--ql',paths_file,'-t',32,'-o',os.path.join(RESULTS_DIR, '{GENUS}/fastani_out')]
