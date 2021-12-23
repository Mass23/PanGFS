import os
import glob
import subprocess

for i in range(0,len(snakemake.input)):
    genomes = glob.glob(snakemake.input[i] + '/*.fasta')
    for genome in genomes:
        name = genome.split('/')[-1].replace('.fasta','')
        args = ['prokka','--cpus',str(snakemake.params.CPU),'--outdir',snakemake.output[i]+'/'+name,'--prefix',name,genome]
        subprocess.call(' '.join(args), shell = True)
