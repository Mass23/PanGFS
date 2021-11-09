# mantis workflow
# Runs MANTIS on the MAGs for metabolic assessment

import os
import glob
import pandas as pd

localrules: bin_collect_mantis, bin_link_mantis, mantis_config, mantis_reformat_consensus, bin_folder_sample_mantis, bin_folder_sample

###########################
# default

rule annotate_genomes:
    input:
        os.path.join(RESULTS_DIR, "logs/mantis.done")
    output:
        touch("status/annotate_genomes.done")


########################################
# checkpoint rules for collecting bins #
########################################
checkpoint bin_collect_mantis:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Binning/selected_DASTool_bins")
    output:
        directory(os.path.join(RESULTS_DIR, "mantis_links/{sample}_bins"))
    shell:
        "ln -vs {input} {output}"

rule bin_link_mantis:
    input:
        os.path.join(RESULTS_DIR, "mantis_links/{sample}_bins/{i}.contigs.fa"),
    output:
        os.path.join(RESULTS_DIR, "mantis_bins/{sample}__{i}.contigs.fa")
    wildcard_constraints:
        sample="|".join(SAMPLES)
        #i="(\w\.)+"
    shell:
        "ln -vs {input} {output}"


####################
# rules for MANTIS #
####################
# MAG list - creating a list of all the mags
# Prokka 
rule prokka:
    input:
        os.path.join(RESULTS_DIR, "mantis_bins/{sample}__{i}.contigs.fa")
    output:
        FAA=os.path.join(RESULTS_DIR, "prokka/{sample}__{i}.faa"),
        FFN=os.path.join(RESULTS_DIR, "prokka/{sample}__{i}.ffn"),
        GFF=os.path.join(RESULTS_DIR, "prokka/{sample}__{i}.gff")
    log:
        os.path.join(RESULTS_DIR, "logs/prokka.{sample}__{i}.log")
    threads:
        config['prokka']['threads']
    conda:
        os.path.join(ENV_DIR, "prokka.yaml")
    message:
        "Running Prokka on {wildcards.mag}"
    shell:
        "(date && prokka --outdir $(dirname {output.FAA}) {input} --cpus {threads} --force && date) &> >(tee {log})"

##########
# MANTIS #
rule mantis_metadata:
    input:
        txt=os.path.join(RESULTS_DIR, "data/mag_list.txt"),
        FAA=rules.prokka.output
    output:
        os.path.join(RESULTS_DIR, "data/mantis_metadata.tsv")
    shell:
        "for fname in {input.txt} ; do echo echo \"${{fname}}\"\"    \"$(echo {input.FAA}) ; done > {output}"

# Mantis: create config file from IMP config
rule mantis_config:
    output:
        "mantis.config"
    message:
        "Creating config for MANTIS"
    run:
        with open(output[0], "w") as ofile:
            # default HMMs
            for hmm_name, hmm_path in config["mantis"]["default"].items():
                ofile.write("%s=%s\n" % (hmm_name, hmm_path))
            # custom HMMs
            for hmm_path in config["mantis"]["custom"]:
                ofile.write("custom_hmm=%s\n" % hmm_path)
            # weights
            for weights_name, weights_value in config["mantis"]["weights"].items():
                ofile.write("%s=%f\n" % (weights_name, weights_value))

# Mantis: protein annotation
# NOTE: check installation before use: python submodules/mantis/ check_installation
rule mantis_run:
    input:
        FAA=rules.prokka.output,
        config="mantis.config"
    output:
        os.path.join(RESULTS_DIR, "mantis/{sample}__{i}/consensus_annotation.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}__{i}.analysis_mantis.log")
    threads:
        config['mantis']['cores']
    params:
         # The below acquires the number of cores to use from the snakemake launch command
#        CORES=int(os.environ.get("CORES"))
        cores=config['mantis']['hmmer_threads']
    conda:
        os.path.join(ENV_DIR, "IMP_MANTIS.yaml")
    message:
        "Annotation with MANTIS for {wildcards.mag}"
    shell:
        "(date && python {config[mantis][path]}/ run_mantis -t {input.FAA} --output_folder $(dirname {output}) --mantis_config {input.config} --hmmer_threads {params.cores} --cores {threads} --memory {config[mantis][single_mem]} --kegg_matrix && date) &> >(tee {log})"

# Mantis: reformat consensus output (to be imported in Python/R)
rule mantis_reformat_consensus:
    input:
        rules.mantis_run.output
    output:
        os.path.join(RESULTS_DIR, "mantis/{sample}__{i}/consensus_annotation.reformatted.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}__{i}.reformat.log")
    message:
        "Reformatting MANTIS output for {wildcards.mag}"
    shell:
        # merge annotations after "|", remove "|" separator
        "(date && paste <(cut -d '|' -f1 {input} | sed 's/\\t$//') <(cut -d '|' -f2 {input} | sed 's/^\\t//;s/\\t/;/g') > {output} && date) &> >(tee {log})"

def bins_mantis(wildcards):
    checkpoint_output = checkpoints.bin_collect_mantis.get(**wildcards).output[0]
    return expand(os.path.join(RESULTS_DIR, "mantis/{sample}__{i}/consensus_annotation.reformatted.tsv"),
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.contigs.fa")).i
        )

rule bin_folder_sample_mantis:
    input:
        bins_mantis
    wildcard_constraints:
        sample="|".join(SAMPLES)
    output:
        touch(os.path.join(RESULTS_DIR, "logs/{sample}_mantis.done"))

rule bin_folder_mantis:
    input:
        expand(os.path.join(RESULTS_DIR, "logs/{sample}_mantis.done"), sample=SAMPLES)
    output:
        touch(os.path.join(RESULTS_DIR, "logs/mantis.done"))

