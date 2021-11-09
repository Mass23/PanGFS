#!/usr/bin/env python
import os
import subprocess


for line in open(snakemake.input[0], 'r').read().split('\n'):
    if line.startswith('GCA'):
        args = ['ncbi-genome-download','-o',snakemake.output[0],'--formats','fasta','--flat-output','-s','genbank','-A',line,'bacteria']     
        subprocess.call(' '.join(args), shell = True)
    elif line.startswith('GCF'):
        args = ['ncbi-genome-download','-o',snakemake.output[0],'--formats','fasta','--flat-output','-s','refseq','-A',line,'bacteria']
        subprocess.call(' '.join(args), shell = True)
    else:
        print(str(line) + ' is not a refseq or genbank accession!')
