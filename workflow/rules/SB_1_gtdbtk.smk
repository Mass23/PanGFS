# mags workflow
# taxonomy of the mags, and then mags selection 

localrules:

###########################
# default

rule gtdbtk:
    input:
        os.path.join(RESULTS_DIR, "gtdbtk_output"),
        os.path.join(DATA_DIR, 'status/cleaning.done')
        #'dummy_test.test'
        #expand(os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs"), GENUS=GENUS_LIST)
        #expand(os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs/{MAG}_clean.fasta"), GENUS=GENUS_LIST)
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
    output:
#        directory(expand(os.path.join(RESULTS_DIR, "{GENUS}/MAGs"), GENUS=GENUS_LIST))
        target=directory(os.path.join(RESULTS_DIR, "{GENUS}/MAGs"))
    run:
        for i in range(0,len(input)):
            os.mkdir(output[i])
            files_to_move = [i for i in open(input[i], 'r').read().split('\n') if i != '']
            for file_to_move in files_to_move:
                contigs_length = [len(rec.seq) for rec in SeqIO.parse(os.path.join(MAGS_DIR, file_to_move) + '.fasta','fasta')]
                if sum(contigs_length) > 100000:
                    os.system('cp -v ' + os.path.join(MAGS_DIR, file_to_move) + '.fasta ' + output[i])

#def aggregate_input(wildcards):
    #checkpoint_output = checkpoints.copy_target_mags.get(**wildcards).output[0]
   # return expand(os.path.join(RESULTS_DIR, '{{GENUS}}/MAGs/{MAG}.fasta'), MAG=glob_wildcards(os.path.join(checkpoint_output, '{MAG}.fasta')).MAG, )
    #files = expand(os.path.join(RESULTS_DIR, '{{GENUS}}/MAGs/{MAG}.fasta'), GENUS=GENUS_LIST, MAG=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)
    #files = [file for file in files if os.path.isfile(file)]
    #return expand(os.path.join(RESULTS_DIR, '{{GENUS}}/MAGs/{MAG}.fasta'), MAG=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)
#    return [file for file in glob.glob(os.path.join(RESULTS_DIR, '*/MAGs/*.fasta'))]

#rule dummy_rule:
#    input:
#        aggregate_input
#    output:
#        'dummy_test.test'
#    shell:
#        'echo {input} > {output}'

rule mag_purify:
    input:
        #aggregate_input
        os.path.join(RESULTS_DIR, "{GENUS}/MAGs/{MAG}.fasta")
    output:
        os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs/{MAG}_clean.fasta")
        #[f.replace('.fasta','_magpurify') for f in glob.glob(os.path.join(RESULTS_DIR, '*/MAGs/*.fasta')],
        #[f.replace('MAGs','cleaned_MAGs').replace('.fasta','_clean.fasta') for f in glob.glob(os.path.join(RESULTS_DIR, '*/MAGs/*.fasta')]
        #directory(os.path.join(RESULTS_DIR, "{GENUS}/MAGs/{MAG}_magpurify")),
        #os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs/{MAG}_clean.fasta")
    wildcard_constraints:
        #MAG=""
        GENUS="|".join(GENUS_LIST)
    conda:
        os.path.join(ENV_DIR, "magpurify.yaml")
    shell:
        "magpurify phylo-markers --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify clade-markers --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify tetra-freq {input} $(basename {input})_magpurify && "
        "magpurify gc-content {input} $(basename {input})_magpurify && "
        "magpurify known-contam --db /mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0 {input} $(basename {input})_magpurify && "
        "magpurify clean-bin {input} $(basename {input})_magpurify {output} " # $(echo {input} | sed -e 's|/MAGS/|/cleaned_MAGs/|g' -e 's|.fasta|_clean.fasta|g') "
    #script:
    #    os.path.join(SRC_DIR, "run_magpurify.py"

def aggregate(wildcards):
    checkpoint_output = checkpoints.copy_target_mags.get(**wildcards).output.target
#    return expand(os.path.dirname(checkpoint_output) + '/cleaned_MAGs/{MAG}_clean.fasta', MAG=glob_wildcards(os.path.dirname(checkpoint_output) + '/cleaned_MAGs/{MAG}_clean.fasta').MAG, )
    return expand(os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs/{MAG}_clean.fasta"), 
        GENUS=GENUS_LIST,
        MAG=glob_wildcards(os.path.join(checkpoint_output, "{MAG}.fasta")).MAG)

rule cleaning_done:
    input:
        aggregate
#    wildcard_constraints:
#        GENUS="|".join(GENUS_LIST)
    output:
        touch(os.path.join(DATA_DIR, 'status/cleaning.done'))

