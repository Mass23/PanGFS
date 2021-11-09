# Pipeline for mapping reads to bins and then identify SNPs
#
# Example call: snakemake -s Snakefile --configfile config.yaml --use-conda --cores 1 -rpn

##############################
# MODULES
import os, re
import glob
import pandas as pd
import numpy as np
#from Bio import SeqIO

##############################
# CONFIG
# can be overwritten by using --configfile <path to config> when calling snakemake
configfile:"config/config.yaml"

##############################
# Relevant directories
DATA_DIR = config["work_dir"]
RESULTS_DIR = config["results_dir"]
ENV_DIR = config["env_dir"]
SRC_DIR = config["scripts_dir"]
MAGS_DIR = config["mags_dir"]
FASTQ_DIR = config["fastq_dir"]

##############################
# Input
GENOMES = "accessions_list.txt"
GENUS_LIST = ['g__Polaromonas','g__Hymenobacter','g__Flavobacterium']