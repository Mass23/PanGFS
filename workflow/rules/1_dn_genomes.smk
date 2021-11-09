# ncbi-genome-download workflow
# Download ref genomes as listed in the accession list

import os
import glob
import pandas as pd

localrules:

###########################
# default 
rule dn_genomes:
    input:
        expand(os.path.join(RESULTS_DIR,'Genomes/{GENUS}'), GENUS=GENUS_LIST)
    output:
        touch("status/dn_genomes.done")

# end it by writing ehe file "genomes_list.txt" that contains the ones that were collected

########################################
# rules for ncbi-download-genomes #
########################################

out_list = []
for key in ACCESSIONS_DICT.keys():
    accessions = ACCESSIONS_DICT[key]
    for accession in accessions:
        out_list.append(os.path.join(RESULTS_DIR,'Genomes',key,accession))

rule download_genomes:
    input:
        GENUS_LIST
    output:
        out_list
    conda:
        os.path.join(ENV_DIR, "ncbi-g-d.yaml")
    script:
        os.path.join(SRC_DIR, "dn_genomes_script.py")

rule create_mags_dir_file:
    input:
        os.path.join(RESULTS_DIR,"Genomes/")
    output:
        os.path.join(DATA_DIR, "genomes_list.txt")
    shell:
        "basename -s '.fna.gz' $(ls {input}/*.fna.gz) > {output}"