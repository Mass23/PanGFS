# Pipeline for MetaPoly
#
# Example call: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores 1 -rpn
# Example call (esb-compute-01 or litcrit): CORES=24 snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${CONDA_PREFIX}/pipeline --cores $CORES -rpn

##############################
# MODULES
import os, re
import glob
import pandas as pd


##############################
# CONFIG
# can be overwritten by using --configfile <path to config> when calling snakemake
# configfile: "config/config.yaml"

include:
    "rules/init.smk"

##############################
# TARGETS & RULES

# List of (main) targets to be created
TARGETS = []
STEPS = ['dn_genomes']#, 'gtdbtk', 'phylo_analysis', 'pan_genome']

# Download genomes
if "dn_genomes" in STEPS:
    include:
        "rules/1_dn_genomes.smk"
    TARGETS += [
        "status/dn_genomes.done"
    ]

# Taxonomy
if "gtdbtk" in STEPS:
    include:
        "rules/2_gtdbtk.smk"
    TARGETS += [
        "status/gtdbtk.done"
    ]

# Phylogeny
if "phylo_analysis" in STEPS:
    include:
        "rules/3_phylo_analysis.smk"
    TARGETS += [
        "status/phylo_analysis.done"
    ]

# PanGenome
if "pan_genome" in STEPS:
    include:
        "rules/4_pan_genome.smk"
    TARGETS += [
        "status/pan_genome.done"
    ]

# No targets
if len(TARGETS) == 0:
    raise Exception('You are not serious. Nothing to be done? Really?')

rule all:
    input:
        TARGETS