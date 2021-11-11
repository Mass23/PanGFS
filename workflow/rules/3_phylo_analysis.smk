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
        expand(os.path.join(RESULTS_DIR, '{GENUS}/gtotree_output/'), GENUS=GENUS_LIST)
    output:
        touch("status/phylo_analysis.done")

rule run_gtotree:
    input:
        GENUS=GENUS_LIST
    output:
        directory(expand(os.path.join(RESULTS_DIR, '{GENUS}/gtotree_output'), GENUS=GENUS_LIST))
    conda:
        os.path.join(ENV_DIR, "gtotree.yaml")
    script:
        os.path.join(SRC_DIR, "run_gtotree.py")
