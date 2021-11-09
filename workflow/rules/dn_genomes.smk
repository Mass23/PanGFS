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
        os.path.join(RESULTS_DIR,"Genomes")
        os.path.join(RESULTS_DIR, "genomes_list.txt")
    output:
        touch("status/dn_genomes.txt")

# end it by writing the file "genomes_list.txt" that contains the ones that were collected

########################################
# rules for ncbi-download-genomes #
########################################

rule download_genomes:
    input:
        GENOMES
    output:
        directory(os.path.join(RESULTS_DIR,"Genomes")),
        os.path.join(RESULTS_DIR,"Genomes/download_genomes.done")
    threads:
        1
    conda:
        os.path.join(ENV_DIR, "ncbi-g-d.yaml")
    run:
        for line in open(input, 'r').readline():
            if line.startswith('GCA'):
                args = ['ncbi-genome-download','-o',output[0],'--formats','fasta','--flat-output','-s','genbank','-A',line,'all']
                subprocess.call(' '.join(args), shell = True)
            elif input.startswith('GCF'):
                args = ['ncbi-genome-download','-o',output[0],'--formats','fasta','--flat-output','-s','refseq','-A',line,'all']
                subprocess.call(' '.join(args), shell = True)
            else:
                print(str(input) + ' is not a refseq or genbank accession!')
        os.system('touch ' + str(output[1]))

rule create_mags_dir:
    input:
        os.path.join(RESULTS_DIR,"Genomes/download_genomes.done")
    output:
        os.path.join(RESULTS_DIR, "genomes_list.txt")
    shell:
        "$(basename -s '.fasta' $(ls $(dirname {input}))/*.fasta) > genomes_list.txt"

