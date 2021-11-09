# ncbi-genome-download workflow
# Download ref genomes as listed in the accession list

import os
import glob
import pandas as pd

localrules:

###########################
# default 

###########################
rule var_calling:
    input:
        ref=os.path.join(BIN_DIR, "{mag}.fa"),
        bam=rules.merge_bam.output
    output:
        gz=os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/{mag}_bam2vcf.log")
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    params:
        calls=temp(os.path.join(RESULTS_DIR, "vcf/{mag}_calls.bcf")),
        filtered=temp(os.path.join(RESULTS_DIR, "vcf/{mag}_filtered.bcf"))
    message:
        "Calling SNPs on {wildcards.mag}"
    shell:
        "(date && mkdir -p $(dirname {output.gz}) && bcftools mpileup -a FORMAT/AD,FORMAT/DP --threads {threads} --max-depth 10000 -f {input.ref} {input.bam} | bcftools call -mv -Ob -o {params.calls} && "
        "bcftools view -i '%QUAL>=30' {params.calls} > {params.filtered} && "
        "bgzip -f {params.filtered} > {output.gz} && date) &> {log}"

