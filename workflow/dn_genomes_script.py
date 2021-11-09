#!/usr/bin/env python
import os
import subprocess

in_genus = snakemake.input[0]

for line in open(os.path.join(DATA_DIR,in_genus) , 'r').read().split('\n'):
    if line.startswith('GCA'):
        args = ['ncbi-genome-download','-o',os.path.join(RESULTS_DIR, 'Genomes', in_genus),'--formats','fasta','--flat-output','-s','genbank','-A',line,'-p',str(snakemake.threads),'bacteria']
        print(' '.join(args))
        subprocess.call(' '.join(args), shell = True)
    elif line.startswith('GCF'):
        args = ['ncbi-genome-download','-o',os.path.join(RESULTS_DIR, 'Genomes', in_genus),'--formats','fasta','--flat-output','-s','refseq','-A',line,'-p',str(snakemake.threads),'bacteria']
        subprocess.call(' '.join(args), shell = True)
    else:
        print(str(line) + ' is not a refseq or genbank accession!')
