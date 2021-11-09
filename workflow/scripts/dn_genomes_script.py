#!/usr/bin/env python
import os
import subprocess
for acc_file in snakemake.input[0]:
    genus=acc_file.split('/')[-2]
    for line in open(acc_file, 'r').read().split('\n'):
        if line.startswith('GCA'):
            args = ['ncbi-genome-download','-o',os.path.join(snakemake.RESULTS_DIR, 'Genomes', genus),'--formats','fasta','--flat-output','-s','genbank','-A',line,'bacteria']     
            subprocess.call(' '.join(args), shell = True)
        elif line.startswith('GCF'):
            args = ['ncbi-genome-download','-o',os.path.join(snakemake.RESULTS_DIR, 'Genomes', genus),'--formats','fasta','--flat-output','-s','refseq','-A',line,'bacteria']
            subprocess.call(' '.join(args), shell = True)
        else:
            print(str(line) + ' is not a refseq or genbank accession!')
