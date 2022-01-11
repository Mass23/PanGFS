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
        'status/dn_genomes.done',
        expand(os.path.join(RESULTS_DIR, '{GENUS}'), GENUS=GENUS_LIST),
        #expand(os.path.join(RESULTS_DIR, '{GENUS}/paths_list_filtered.txt'), GENUS=GENUS_LIST)
        expand(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'), GENUS=GENUS_LIST),
        #expand(os.path.join(RESULTS_DIR, '{GENUS}/Mappings'), GENUS=GENUS_LIST),
        #expand(os.path.join(RESULTS_DIR, '{GENUS}/pangenome_annotations_list.txt'), GENUS=GENUS_LIST),
        expand(os.path.join(RESULTS_DIR, '{GENUS}/mmseqs2_all_seqs.fasta'), GENUS=GENUS_LIST)
        #os.path.join(RESULTS_DIR, "logs/mmseqs.done")
    output:
        touch("status/pan_genome.done")

#rule dn_outgroup:
#    input:
#        os.path.join(DATA_DIR, 'status/phylo_analysis.done')
#    output:
#        os.path.join(RESULTS_DIR, '{GENUS}/Genomes/outgroup.fasta')

rule prokka_genomes:
    input:
        expand(os.path.join(RESULTS_DIR, '{GENUS}'), GENUS=GENUS_LIST)
    output:
        directory(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'))
    params:
        CPU=12
    conda:
        os.path.join(ENV_DIR, "prokka.yaml")
    script:
        os.path.join(SRC_DIR, "run_prokka_genomes.py")

#rule prokka_mags:
#    input:
#        expand(os.path.join(RESULTS_DIR, '{GENUS}/cleaned_MAGs'), GENUS=GENUS_LIST)
#    output:
#        directory(expand(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'), GENUS=GENUS_LIST))
#    params:
#        CPU=32
#    conda:
#        os.path.join(ENV_DIR, "prokka.yaml")
#    script:
#        os.path.join(SRC_DIR, "run_prokka_mags.py")

#rule list_genome_annotations:
#    input:
#        expand(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'), GENUS=GENUS_LIST)
#    output:
#        expand(os.path.join(RESULTS_DIR, '{GENUS}/pangenome_annotations_list.txt'), GENUS=GENUS_LIST)
#    run:
#        import os
#        for i in range(0,len(input)):
#            os.system('ls ' + input[i] + '/*/GC*.faa > ' + output[i])

#rule run_mmseqs2_cluster:
#    input:
#        expand(os.path.join(RESULTS_DIR, '{GENUS}/Pangenome'), GENUS=GENUS_LIST)
#    output:
#        expand(os.path.join(RESULTS_DIR, '{GENUS}/mmseqs2_all_seqs.fasta'), GENUS=GENUS_LIST)
#    conda:
#        os.path.join(ENV_DIR, "mmseqs2.yaml")
#    script:
#        os.path.join(SRC_DIR, 'run_mmseqs2.py')




