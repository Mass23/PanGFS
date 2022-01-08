import os
import glob
import subprocess
def process_genus(in_folder, out_folder):
    print(in_folder)
    print(out_folder)
    raw_mags = glob.glob(in_folder + '/*.fasta')
    clean_mags = [i.split('/')[-1].replace('.fasta','_clean.fasta') for i in raw_mags]
    for mag in range(0,len(raw_mags)):
        mag_in = os.path.join(in_folder, raw_mags[mag])
        mag_puri_out = os.path.join(in_folder, raw_mags[mag].replace('.fasta','_magpurify'))
        mag_out = out_folder + '/' + clean_mags[mag]
        args1 = ['magpurify', 'phylo-markers',
        '--db', '/mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0', mag_in, mag_puri_out]
        subprocess.call(' '.join(args1), shell = True)
        args2 = ['magpurify', 'clade-markers',
        '--db', '/mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0', mag_in, mag_puri_out]
        subprocess.call(' '.join(args2), shell = True)
        args3 = ['magpurify', 'tetra-freq', mag_in, mag_puri_out]
        subprocess.call(' '.join(args3), shell = True)
        args4 = ['magpurify', 'gc-content', mag_in, mag_puri_out]
        subprocess.call(' '.join(args4), shell = True)
        args5 = ['magpurify', 'known-contam',
        '--db', '/mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0', mag_in, mag_puri_out]
        subprocess.call(' '.join(args5), shell = True)
        args_clean = ['magpurify', 'clean-bin', mag_in, mag_puri_out, mag_out]
        subprocess.call(' '.join(args_clean), shell = True)

#for i in range(0,len(snakemake.input)):
os.makedirs(str(snakemake.output))
process_genus(str(snakemake.input), str(snakemake.output))
