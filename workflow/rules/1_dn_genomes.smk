# ncbi-genome-download workflow
# Download ref genomes as listed in the accession list

import os
import glob
import pandas as pd

localrules:

###########################
rule dn_genomes:
    input:
        expand(RESULTS_DIR + '/Genomes/{GENUS}/',GENUS=GENUS_LIST),
        expand(os.path.join(DATA_DIR, "{GENUS}/genomes_list.txt"), GENUS=GENUS_LIST)
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

rule create_genomes_list:
    input:
        expand(os.path.join(RESULTS_DIR,"Genomes/{GENUS}"), GENUS=GENUS_LIST)
    output:
        expand(os.path.join(DATA_DIR, "{GENUS}/genomes_list.txt"), GENUS=GENUS_LIST)
    run:
        import os
        for i in range(0,len(snakemake.input)):
            os.system('gunzip ' + snakemake.input[i] + '/*.fna.gz')
            os.system('ls ' + snakemake.input[i] + '/*.fna > ' + snakemake.output[i])
