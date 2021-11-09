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
GENUS_LIST = ['g__Polaromonas','g__Hymenobacter','g__Flavobacterium','g__Methylotenera','g__Novosphingobium']
ACCESSIONS_DICT = {'g__Polaromonas': open(os.path.join(DATA_DIR,'accessions/Polaromonas/Polaromonas_accessions'),'r').read().split('\n'),
                   'g__Hymenobacter': open(os.path.join(DATA_DIR,'accessions/Hymenobacter/Hymenobacter_accessions'),'r').read().split('\n'),
                   'g__Flavobacterium': open(os.path.join(DATA_DIR,'accessions/Flavobacterium/Flavobacterium_accessions'),'r').read().split('\n'),
                   'g__Methylotenera': open(os.path.join(DATA_DIR,'accessions/Methylotenera/Methylotenera_accessions'),'r').read().split('\n'),
                   'g__Novosphingobium': open(os.path.join(DATA_DIR,'accessions/Novosphingobium/Novosphingobium_accessions'),'r').read().split('\n')}