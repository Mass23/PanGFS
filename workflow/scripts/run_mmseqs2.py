import os
import subprocess

for i in range(0,len(snakemake.input)):
    args = ['mmseqs', 'easy-cluster', snakemake.input[i] + '/*/*.faa', snakemake.output[i].replace('_all_seqs.fasta',''), 'tmp_cluster', '--min-seq-id', '0.3', '-c', '0.8', '--cov-mode', '0', '--threads', '32']
    subprocess.call(' '.join(args), shell=True)
