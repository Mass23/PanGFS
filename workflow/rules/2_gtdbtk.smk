# mags workflow
# taxonomy of the mags, and then mags selection 

localrules:

###########################
# default

rule gtdbtk:
    input:
        os.path.join(RESULTS_DIR, "gtdbtk_output/"),
        expand(os.path.join(DATA_DIR, "{GENUS}/mags_list.txt"),GENUS=GENUS_LIST),
        expand(os.path.join(RESULTS_DIR, "{GENUS}/MAGs/"), GENUS=GENUS_LIST),
        expand(os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs/"), GENUS=GENUS_LIST)
    output:
        touch("status/gtdbtk.done")


################################################
# rule for collecting mags with taxo. matching #
################################################

rule run_gtdbtk:
    input:
        os.path.join(MAGS_DIR)
    output:
        directory(os.path.join(RESULTS_DIR, "gtdbtk_output/"))
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
        GTDBTK_DIR=os.path.join(RESULTS_DIR, "gtdbtk_output/"),
        GENUS=GENUS_LIST
    output:
        expand(os.path.join(DATA_DIR, "{GENUS}/mags_list.txt"),GENUS=GENUS_LIST)
    run:
        for i in range(0,len(GENUS_LIST)):
            tax_string = input.GENUS[i]
            gtdbtk_file = pd.read_csv(input.GTDBTK_DIR + '/gtdbtk.bac120.summary.tsv', sep='\t')
            gtdbtk_sub = gtdbtk_file[gtdbtk_file['classification'].str.contains(tax_string)].user_genome
            gtdbtk_sub.to_csv(output[i], header =  False, sep='\t', index=False)

rule copy_target_mags:
    input:
        expand(os.path.join(DATA_DIR, "{GENUS}/mags_list.txt"),GENUS=GENUS_LIST)
    output:
        directory(expand(os.path.join(RESULTS_DIR, "{GENUS}/MAGs"), GENUS=GENUS_LIST))
    run:
        for i in range(0,len(input)):
            os.mkdir(output[i])
            files_to_move = open(input[i], 'r').read().split('\n')
            for file_to_move in files_to_move:
                contigs_length = [len(rec.seq) for rec in SeqIO.parse(file_to_move,'fasta')]
                if sum(contigs_length) > 100000:
                    os.system('cp -v ' + os.path.join(MAGS_DIR, file_to_move) + '.fasta ' + output[i])

rule mag_purify:
    input:
        expand(os.path.join(RESULTS_DIR, "{GENUS}/MAGs/"), GENUS=GENUS_LIST)
    output:
        directory(expand(os.path.join(RESULTS_DIR, "{GENUS}/cleaned_MAGs"), GENUS=GENUS_LIST))
    run:
        os.system('conda activate magpurify')
        os.system('export MAGPURIFYDB=/mnt/esb-storage-01/NOMIS/databases/MAGpurify-db-v1.0')
        for i in range(0,len(input)):
            in_folder = input[i]
            out_folder = output[i]
            raw_mags = glob.glob(in_folder)
            clean_mags = [i.replace('.fasta','_clean.fasta') for i in raw_mags]
            for mag in range(0,len(raw_mags)):
                mag_in = os.path.join(in_folder, raw_mags[i])
                mag_puri_out = os.path.join(in_folder, raw_mags[i].replace('.fasta','_magpurify/'))
                mag_out = os.path.join(out_folder, clean_mags[i])
                args1 = ['magpurify', 'phylo-markers', mag_in, mag_puri_out]
                subprocess.call(' '.join(args1), shell = True)
                args2 = ['magpurify', 'clade-markers', mag_in, mag_puri_out]
                subprocess.call(' '.join(args2), shell = True)
                args3 = ['magpurify', 'tetra-freq',    mag_in, mag_puri_out]
                subprocess.call(' '.join(args3), shell = True)
                args4 = ['magpurify', 'gc-content',    mag_in, mag_puri_out]
                subprocess.call(' '.join(args4), shell = True)
                args5 = ['magpurify', 'known-contam',  mag_in, mag_puri_out]
                subprocess.call(' '.join(args5), shell = True)
                args_clean = ['magpurify', 'clean-bin', mag_in, mag_puri_out, mag_out]
                subprocess.call(' '.join(args6), shell = True)


