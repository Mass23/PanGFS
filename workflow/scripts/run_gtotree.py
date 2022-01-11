#!/usr/bin/env python
import os
import subprocess
import glob

def process_genus(genus, hmm, out, res_dir):
    # list paths
    pan_annot_files = os.path.join(res_dir, genus, 'Pangenome/*/*.faa')
    pan_annot_list = os.path.join(res_dir, genus, 'Pangenome/genomes_list.txt')
    os.system('ls ' + pan_annot_files + ' > ' + pan_annot_list)

    # create file containing paths/accessions
    acc_out = os.path.join(res_dir, genus, 'acc_outgroup.txt')
    with open(acc_out, 'w') as f:
        f.write(out)
    # run  gtotree
    gtotree_args = ['GToTree','-o',os.path.join(res_dir,genus,'gtotree_out'),
                    '-A',pan_annot_list,
                    '-a',acc_out,'-N',
                    '-H',hmm,'-j','48','-G','0.3']
    subprocess.call(' '.join(gtotree_args), shell = True)
    #iqtree_args = ['iqtree','-s',os.path.join(res_dir,genus,'gtotree_out/Aligned_SCGs.faa'),
    #                '-spp',os.path.join(res_dir, genus,'gtotree_out/run_files/Partitions.txt'),
    #                '-m','MFP','-nt','48','-bb','1000','-pre',os.path.join(res_dir, genus,'iqtree_out')]
    #subprocess.call(' '.join(iqtree_args), shell = True)

for i in range(0,len(snakemake.params.GENUS)):
    process_genus(snakemake.params.GENUS[i], snakemake.params.HMM[i], snakemake.params.OUT[i], snakemake.params.RES_DIR)
