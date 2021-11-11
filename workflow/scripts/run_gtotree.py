#!/usr/bin/env python
import os
import subprocess
import glob

def process_genus(genus, hmm, out, res_dir):
    # list paths
    genomes_dir = os.path.join(res_dir, 'Genomes', genus)
    genomes_paths = [os.path.join(res_dir, 'Genomes', genus,gen_path) for gen_path in os.listdir(genomes_dir) if gen_path.endswith('.fna') or gen_path.endswith('.fasta')]
    mags_dir = os.path.join(res_dir, 'MAGs', genus)
    mags_paths = [os.path.join(res_dir, 'MAGs', genus, mag_path) for mag_path in os.listdir(mags_dir) if mag_path.endswith('.fna') or mag_path.endswith('.fasta')]
    merged_paths = genomes_paths + mags_paths
    # create file containing paths/accessions
    with open(os.path.join(res_dir, genus, 'paths_list.txt'), 'w') as f:
        f.write('\n'.join(merged_paths))
    with open(os.path.join(res_dir, genus, 'acc_outgroup.txt'), 'w') as f:
        f.write(out)
    # run  gtotree
    gtotree_args = ['GToTree','-o',os.path.join(res_dir,genus,'gtotree_out'),
                    '-f',os.path.join(res_dir, genus, 'paths_list.txt'),
                    '-a',os.path.join(res_dir, genus, 'acc_outgroup.txt'),
                    '-H',hmm,'-j','32']
    subprocess.call(' '.join(gtotree_args), shell = True)

for i in range(0,len(snakemake.input.GENUS)):
    process_genus(snakemake.input.GENUS[i], snakemake.params.HMM[i], snakemake.params.OUT[i], snakemake.params.RES_DIR)