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
    run:
        os.system('conda activate gtotree')
        for i in range(0,len(GENUS_LIST)):
            # list paths
            genomes_dir = os.path.join(RESULTS_DIR, 'Genomes', input.GENUS[i])
            genomes_paths = [os.path.join(RESULTS_DIR, 'Genomes', input.GENUS[i],gen_path) for gen_path in genomes_dir]
            mags_dir = os.path.join(RESULTS_DIR, 'MAGs', input.GENUS[i])
            mags_paths = [os.path.join(RESULTS_DIR, 'MAGs', input.GENUS[i],mag_path) for mag_path in mags_dir]
            merged_paths = genomes_paths + mags_paths
            # create file containing paths/accessions
            with open(os.path.join(RESULTS_DIR, input.GENUS[i], 'paths_list.txt')) as f:
                f.write('\n'.join(merged_paths))
            with open(os.path.join(RESULTS_DIR, input.GENUS[i], 'acc_outgroup.txt')) as f:
                f.write(OUTGROUP_DICT[input.GENUS[i]])
            # run  gtotree
            gtotree_args = ['GToTree','-o',os.path.join(RESULTS_DIR,input.GENUS[i],'gtotree_out'),
                            '-f',os.path.join(RESULTS_DIR, input.GENUS[i], 'paths_list.txt'),
                            '-a',os.path.join(RESULTS_DIR, input.GENUS[i], 'acc_outgroup.txt'),
                            '-H',HMM_DICT[input.GENUS[i]],'-j','32']
            subprocess.call(' '.join(gtotree_args), shell = True)
