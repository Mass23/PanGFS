import os
import subprocess

for i in range(0,len(snakemake.input)):
    args = ['mmseqs', 'easy-cluster', snakemake.input[i] + '/*/*.faa', snakemake.output[i].replace('_all_seqs.fasta',''), 'tmp_cluster', '--min-seq-id', '0.7', '-c', '0.7', '--cov-mode', '1', '--threads', '32']
    subprocess.call(' '.join(args), shell=True)
