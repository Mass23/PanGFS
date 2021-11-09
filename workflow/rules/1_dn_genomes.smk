# ncbi-genome-download workflow
# Download ref genomes as listed in the accession list

import os
import glob
import pandas as pd

localrules:

###########################
# default 
#out_list = []
#for key in ACCESSIONS_DICT.keys():
#    acc_file =  ACCESSIONS_DICT[key]
#    acc_list = open(acc_file, 'r').read().split('\n')
#    for accession in acc_list:
#        out_list.append(os.path.join(RESULTS_DIR,'Genomes',key,accession+'_{wildcards.name}.fna.gz'))

rule dn_genomes:
    input:
        expand(RESULTS_DIR + '/Genomes/{GENUS}/',GENUS=GENUS_LIST)
    output:
        touch("status/dn_genomes.done")

# end it by writing ehe file "genomes_list.txt" that contains the ones that were collected

########################################
# rules for ncbi-download-genomes #
########################################

rule download_genomes:
    input:
        expand(DATA_DIR + '/accessions/{GENUS}/{GENUS}_accessions',GENUS=GENUS_LIST)
    output:
        directory(expand(RESULTS_DIR + '/Genomes/{GENUS}/',GENUS=GENUS_LIST))
    conda:
        os.path.join(ENV_DIR, "ncbi-g-d.yaml")
    script:
        os.path.join(SRC_DIR, "dn_genomes_script.py")
