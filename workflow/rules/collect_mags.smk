# bins workflow
# Collects all the bins from each sample in single folder 

localrules: bin_collect, bin_link, bin_folder_sample, bin_folder

###########################
# default

rule bins:
    input:
        os.path.join(RESULTS_DIR, "data/mag_list.txt")
    output:
        touch("status/bins.done")


########################################
# checkpoint rules for collecting bins #
########################################
checkpoint bin_collect:
    input:
        os.path.join(DATA_DIR, "{sample}/run1/Binning/selected_DASTool_bins")
    output:
        directory(os.path.join(RESULTS_DIR, "links/{sample}_bins"))
    shell:
        "ln -vs {input} {output}"

rule bin_link:
    input:
        os.path.join(RESULTS_DIR, "links/{sample}_bins/{i}.contigs.fa"),
    output:
        os.path.join(RESULTS_DIR, "bins/{sample}__{i}.contigs.fa")
    wildcard_constraints:
        sample="|".join(SAMPLES)
        #i="(\w\.)+"
    shell:
        "ln -vs {input} {output}"

def bins(wildcards):
    checkpoint_output = checkpoints.bin_collect.get(**wildcards).output[0]
    return expand(os.path.join(RESULTS_DIR, "bins/{{sample}}__{i}.contigs.fa"),
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.contigs.fa")).i)

rule bin_folder_sample:
    input:
        bins
    wildcard_constraints:
        sample="|".join(SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "logs/{sample}_bin_collection.done"),
    shell:
        "for fname in {input} ; do echo $(basename -s \".contigs.fa\" \"${{fname}}\") ; done > {output}"
        

rule bin_folder:
    input:
        expand(os.path.join(RESULTS_DIR, "logs/{sample}_bin_collection.done"), sample=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "data/mag_list.txt"),
        touch(os.path.join(RESULTS_DIR, "bins/bin_collection.done"))
    shell:
        "cat {input} > {output[0]}"
