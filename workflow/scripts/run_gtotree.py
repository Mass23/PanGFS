#!/usr/bin/env python
import os
import subprocess
import glob

def process_genus(genus, hmm, out, res_dir):
    # list paths
    genomes_dir = os.path.join(res_dir, genus, 'Genomes')
    genomes_paths = [os.path.join(res_dir, genus, 'Genomes',gen_path) for gen_path in os.listdir(genomes_dir) if gen_path.endswith('.fna') or gen_path.endswith('.fasta')]
    mags_dir = os.path.join(res_dir, genus, 'cleaned_MAGs')
    mags_paths = [os.path.join(res_dir, genus, 'cleaned_MAGs', mag_path) for mag_path in os.listdir(mags_dir) if mag_path.endswith('.fna') or mag_path.endswith('.fasta')]
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
                    '-H',hmm,'-j','32','-G','0.5','-N']
    subprocess.call(' '.join(gtotree_args), shell = True)
    iqtree_args = ['iqtree','-s',os.path.join(res_dir,genus,'gtotree_out/Aligned_SCGs.faa'),
                    '-spp',os.path.join(res_dir, genus,'gtotree_out/run_files/Partitions.txt'),
                    '-m','MFP','-nt','32','-bb','1000','-pre',os.path.join(res_dir, genus,'iqtree_')]
    subprocess.call(' '.join(iqtree_args), shell = True)

for i in range(0,len(snakemake.params.GENUS)):
    process_genus(snakemake.params.GENUS[i], snakemake.params.HMM[i], snakemake.params.OUT[i], snakemake.params.RES_DIR)