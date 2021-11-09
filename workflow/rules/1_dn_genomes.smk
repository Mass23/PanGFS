# ncbi-genome-download workflow
# Download ref genomes as listed in the accession list

import os
import glob
import pandas as pd

localrules:

###########################
# default 
out_list = []
for key in ACCESSIONS_DICT.keys():
    acc_file =  ACCESSIONS_DICT[key]
    acc_list = open(acc_file, 'r').read().split('\n')
    for accession in acc_list:
        out_list.append(os.path.join(RESULTS_DIR,'Genomes',key,accession+'_{wildcards.name}.fna'))

rule dn_genomes:
    input:
        out_list
    output:
        touch("status/dn_genomes.done")

# end it by writing ehe file "genomes_list.txt" that contains the ones that were collected

########################################
# rules for ncbi-download-genomes #
########################################

rule download_genomes:
    input:
        [ACCESSIONS_DICT[key] for key in ACCESSIONS_DICT.keys()]
    output:
        out_list
    conda:
        os.path.join(ENV_DIR, "ncbi-g-d.yaml")
    run:
        import os
        import subprocess

        acc_file = input[0]
        genus=acc_file.split('/')[-2]

        for line in open(acc_file, 'r').read().split('\n'):
            if line.startswith('GCA'):
                args = ['ncbi-genome-download','-o',os.path.join(RESULTS_DIR, 'Genomes', genus),'--formats','fasta','--flat-output','-s','genbank','-A',line,'bacteria']
                subprocess.call(' '.join(args), shell = True)
            elif line.startswith('GCF'):
                args = ['ncbi-genome-download','-o',os.path.join(RESULTS_DIR, 'Genomes', genus),'--formats','fasta','--flat-output','-s','refseq','-A',line,'bacteria']
                subprocess.call(' '.join(args), shell = True)
            else:
                print(str(line) + ' is not a refseq or genbank accession!')
