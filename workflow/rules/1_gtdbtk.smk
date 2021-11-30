# mags workflow
# taxonomy of the mags, and then mags selection 

localrules:

###########################
# default

rule gtdbtk:
    input:
        os.path.join(RESULTS_DIR, "gtdbtk_output/")#,
        #expand(os.path.join(DATA_DIR, "accessions/{GENUS}/mags_list.txt"),GENUS=GENUS_LIST),
        #expand(os.path.join(RESULTS_DIR, "{GENUS}/MAGs/"), GENUS=GENUS_LIST),
        #expand(os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs/"), GENUS=GENUS_LIST)
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

rule copy_target_mags:
    input:
        expand(os.path.join(DATA_DIR, "accessions/{GENUS}/mags_list.txt"),GENUS=GENUS_LIST)
    output:
        directory(expand(os.path.join(RESULTS_DIR, "{GENUS}/MAGs"), GENUS=GENUS_LIST))
    run:
        for i in range(0,len(input)):
            os.mkdir(output[i])
            files_to_move = [i for i in open(input[i], 'r').read().split('\n') if i != '']
            for file_to_move in files_to_move:
                contigs_length = [len(rec.seq) for rec in SeqIO.parse(os.path.join(MAGS_DIR, file_to_move) + '.fasta','fasta')]
                if sum(contigs_length) > 100000:
                    os.system('cp -v ' + os.path.join(MAGS_DIR, file_to_move) + '.fasta ' + output[i])

rule mag_purify:
    input:
        expand(os.path.join(RESULTS_DIR, "{GENUS}/MAGs/"), GENUS=GENUS_LIST)
    output:
        directory(expand(os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs"), GENUS=GENUS_LIST))
    conda:
        os.path.join(ENV_DIR, "magpurify.yaml")
    script:
        os.path.join(SRC_DIR, "run_magpurify.py")