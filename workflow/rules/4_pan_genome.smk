# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: bin_collect_mantis, bin_link_mantis, mantis_config, mantis_reformat_consensus, bin_folder_sample_mantis, bin_folder_sample

###########################
# default

rule pan_genome:
    input:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/paths_list_filtered.txt'), GENUS=GENUS_LIST)
        #os.path.join(RESULTS_DIR, "logs/prokka.done"),
        #os.path.join(RESULTS_DIR, "logs/mmseqs.done")
    output:
        touch("status/annotate_genomes.done")


rule filter_paths_list:
    input:
        PATH_F=expand(os.path.join(RESULTS_DIR, '{GENUS}/paths_list.txt'), GENUS=GENUS_LIST),
        TO_RM=expand(os.path.join(RESULTS_DIR, '{GENUS}/gtotree_out/Genomes_removed_for_too_few_hits.tsv'), GENUS=GENUS_LIST)
    output:
        expand(os.path.join(RESULTS_DIR, '{GENUS}/paths_list_filtered.txt'), GENUS=GENUS_LIST)
    run:
        genomes_to_remove = pd.read_csv(TO_RM, sep='\t', header=True)['accession'].tolist()
        with open(input.PATH_F, 'r') as in_file:
            with open(output) as out_file:
                for line in in_file.readline():
                    if line.split('/')[-1] in genomes_to_remove:
                        continue
                    else:
                        out_file.write(line)
