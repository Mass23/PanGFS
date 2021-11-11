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
        GENUS=GENUS_LIST,
        OUTGROUP=OUTGROUP_LIST,
        HMM_gtotree=HMM_GTOTREE
    output:
        expand(directory(os.path.join(RESULTS_DIR, "{GENUS}/gtotree_output")), GENUS=GENUS_LIST)
    conda:
        os.path.join(ENV_DIR, "gtotree.yaml")
    params:
        hmm_gtotree=config['gtotree']['hmm'],
        threads_gtotree=config['gtotree']['threads']
    run:
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
                f.write(input.OUTGROUP[i])
            gtotree_args = ['GToTree','-o',os.path.join(RESULTS_DIR,input.GENUS[i],'gtotree_out'),
                            '-f',os.path.join(RESULTS_DIR, input.GENUS[i], 'paths_list.txt'),
                            '-a',os.path.join(RESULTS_DIR, input.GENUS[i], 'acc_outgroup.txt'),
                            '-H',input.HMM_gtotree[i],'-j','32']
