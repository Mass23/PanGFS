#!/usr/bin/env python
import os
import subprocess

def process_genus(genus_acc, genus_dir):
    for line in open(genus_acc, 'r').read().split('\n'):
        if line.startswith('GCA'):
            args = ['ncbi-genome-download','-o',genus_dir,'--formats','fasta','--flat-output','-s','genbank','--parallel','8','-A',line,'bacteria']     
            subprocess.call(' '.join(args), shell = True)
        elif line.startswith('GCF'):
            args = ['ncbi-genome-download','-o',genus_dir,'--formats','fasta','--flat-output','-s','refseq','--parallel','8','-A',line,'bacteria']
            subprocess.call(' '.join(args), shell = True)
        else:
            print(str(line) + ' is not a refseq or genbank accession!')

for i in range(0,len(snakemake.input)):
    process_genus(snakemake.input[i],snakemake.output[i])