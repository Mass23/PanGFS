# Pipeline for mapping reads to bins and then identify SNPs
#
# Example call: snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --cores 1 -rpn

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

##############################
# Input
GENUS_LIST = ['g__Polaromonas','g__Hymenobacter','g__Flavobacterium','g__Methylotenera','g__Novosphingobium']
ACCESSIONS_DICT = {'g__Polaromonas': os.path.join(DATA_DIR,'accessions/Polaromonas/Polaromonas_accessions'),
                   'g__Hymenobacter': os.path.join(DATA_DIR,'accessions/Hymenobacter/Hymenobacter_accessions'),
                   'g__Flavobacterium': os.path.join(DATA_DIR,'accessions/Flavobacterium/Flavobacterium_accessions'),
                   'g__Methylotenera': os.path.join(DATA_DIR,'accessions/Methylotenera/Methylotenera_accessions'),
                   'g__Novosphingobium': os.path.join(DATA_DIR,'accessions/Novosphingobium/Novosphingobium_accessions')}