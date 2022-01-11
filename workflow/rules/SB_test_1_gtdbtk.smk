# mags workflow
# taxonomy of the mags, and then mags selection 

localrules:

###########################
# default

rule gtdbtk:
    input:
        os.path.join(RESULTS_DIR, "gtdbtk_output"),
        os.path.join(DATA_DIR, 'status/cleaning.done')
    output:
        touch("status/gtdbtk.done")


################################################
# rule for collecting mags with taxo. matching #
################################################

rule run_gtdbtk:
    input:
        os.path.join(MAGS_DIR)
    output:
        directory(os.path.join(RESULTS_DIR, "gtdbtk_output"))
    log:
        os.path.join(RESULTS_DIR, "logs/gtdbtk.log")
    conda:
        os.path.join(ENV_DIR, "gtdbtk_updated.yaml")
    params:
        config["gtdbtk"]["path"]
    threads:
        config["gtdbtk"]["threads"]
    message:
        "Running GTDB toolkit on MAGs"
    shell:
        "export OMP_NUM_THREADS=48 && "
        "export PENBLAS_NUM_THREADS=48 && "
        "export MKL_NUM_THREADS=48 && "
        "export VECLIB_MAXIMUM_THREADS=48 && "
        "export NUMEXPR_NUM_THREADS=48 && "
        "(date && export GTDBTK_DATA_PATH={params} && gtdbtk classify_wf --cpus {threads} -x fasta --genome_dir {input} --out_dir {output} && date) &> {log}"

rule list_target_mags:
    input:
        GTDBTK_DIR=os.path.join(RESULTS_DIR, "gtdbtk_output/")
    params:
        GENUS=GENUS_LIST
    output:
        expand(os.path.join(DATA_DIR, "accessions/{GENUS}/mags_list.txt"),GENUS=GENUS_LIST)
    run:
        for i in range(0,len(GENUS_LIST)):
            tax_string = 'g__' + params.GENUS[i]
            gtdbtk_file = pd.read_csv(input.GTDBTK_DIR + '/gtdbtk.bac120.summary.tsv', sep='\t')
            gtdbtk_sub = gtdbtk_file[gtdbtk_file['classification'].str.contains(tax_string)].user_genome
            gtdbtk_sub.to_csv(output[i], header =  False, sep='\t', index=False)

checkpoint copy_target_mags:
    input:
        expand(os.path.join(DATA_DIR, "accessions/{GENUS}/mags_list.txt"),GENUS=GENUS_LIST)
#        os.path.join(DATA_DIR, "accessions/{GENUS}/mags_list.txt")
    output:
        target=directory(expand(os.path.join(RESULTS_DIR, "{GENUS}/MAGs"), GENUS=GENUS_LIST))
    run:
        for i in range(0,len(input)):
            os.mkdir(output[i])
            files_to_move = [i for i in open(input[i], 'r').read().split('\n') if i != '']
            for file_to_move in files_to_move:
                contigs_length = [len(rec.seq) for rec in SeqIO.parse(os.path.join(MAGS_DIR, file_to_move) + '.fasta','fasta')]
                if sum(contigs_length) > 100000:
                    os.system('cp -v ' + os.path.join(MAGS_DIR, file_to_move) + '.fasta ' + output[i])

rule bin_link:
    input:
        os.path.join(RESULTS_DIR, "{GENUS}/MAGs/{i}.fasta")
    output:
        os.path.join(RESULTS_DIR, "linked/{GENUS}/MAGs/{i}.fasta")
    wildcard_constraints:
        GENUS="|".join(GENUS_LIST)
    shell:
        "ln -vs {input} {output}"

rule mag_purify:
    input:
        os.path.join(RESULTS_DIR, "linked/{GENUS}/MAGs/{i}.fasta")
    output:
        os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs/{i}_clean.fasta")
        #os.path.join(RESULTS_DIR, "magpurify/{GENUS}/cleaned_MAGs/{i}_clean.fasta")
    wildcard_constraints:
        GENUS="|".join(GENUS_LIST)
    conda:
        os.path.join(ENV_DIR, "magpurify.yaml")
    shell:
        "magpurify phylo-markers --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} {input}_magpurify && "
        "magpurify clade-markers --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} {input}_magpurify && "
        "magpurify tetra-freq {input} {input}_magpurify && "
        "magpurify gc-content {input} {input}_magpurify && "
        "magpurify known-contam --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} {input}_magpurify && "
        #"magpurify clean-bin {input} $(basename {input})_magpurify {output} "  # $(echo {input} | sed -e 's|/MAGS/|/cleaned_MAGs/|g' -e 's|.fasta|_clean.fasta|g') "
        "magpurify clean-bin {input} {input}_magpurify $(basename {output}) " 

def aggregate_mags(wildcards):
    checkpoint_output = checkpoints.copy_target_mags.get(**wildcards).output.target
    print(checkpoint_output)
    print(glob_wildcards(os.path.join(checkpoint_output[0], "{i}.fasta")))
    mags=[]
    for genus_folder in checkpoint_output:
        genus=os.path.basename(os.path.dirname(genus_folder))
        mags.extend(
            expand(
                os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs/{i}.fasta"), GENUS=genus, i=glob_wildcards(os.path.join(genus_folder, "{i}.fasta")).i
            )
        )
    return mags

#    return expand(os.path.join(RESULTS_DIR, "magpurify/{GENUS}/cleaned_MAGs/{i}_clean.fasta"), 
#        GENUS=GENUS_LIST,
#        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)

rule cleaning_done:
    input:
        aggregate_mags
#    wildcard_constraints:
#        GENUS="|".join(GENUS_LIST)
    output:
        touch(os.path.join(DATA_DIR, 'status/cleaning.done'))
