# iDSBseq : Identification of DNA double-strand breaks
#
# Author : Christopher Nobles, Ph.D.
# useful modules sys, os, csv

import os
import sys
import re
import yaml
import configparser
from tools.pytools.defs import *

if not config:
    raise SystemExit(
        "No config file specified. Run `sunbeam_init` to generate a "
        "config file, and specify with --configfile")

# Import sampleInfo
if ".csv" in config["sampleInfo"]:
    delim = ","
elif ".tsv" in config["sampleInfo"]:
    delim = "\t"
else:
    raise SystemExit("Sample Info file contain extention '.csv' or '.tsv'.")

# Sample information
sampleInfo = import_sample_info(
    config["Sample_Info"], config["Sample_Name_Column"], delim)

SAMPLES=sampleInfo[config["sample_name_column"]]
TYPES=config["Read_Types"]
READS=config["Genomic_Reads"]

# Working paths
RUN = config["Run_Name"]
ROOT_DIR = config["Install_Directory"]
RUN_DIR = config["Install_Directory"] + "/analysis/" + RUN

# Check for directory paths !!! Not sure if this is going to work just yet.
if not os.path.isdir(ROOT_DIR):
    raise SystemExit("Path to iDSBseq is not found. Check configuration file.")
if not os.path.isdir(RUN_DIR):
    raise SystemExit(
        "Analysis directory is not constructed for the run.",
        "Please run 'snakemake setup_working_dir -c [config.file]'.")

# Architecture Rules
include: "rules/workflow_misc/arch.rules"

# Processing Rules
include: "rules/demultiplex/demulti.rules"
include: "rules/sequence_trim/trim.rules"
include: "rules/consolidate/consol.rules"
if (config["Aligner"] == "BLAT" or config["Aligner"] == "blat"):
    include: "rules/align/align.blat.rules"
    include: "rules/post_align/post_align.blat.rules"
elif (config["Aligner"] == "BWA" or config["Aligner"] == "bwa"):
    raise SystemExit("BWA aligner not supported yet.")
else:
    "Aligner: " + config["Aligner"] + " not supported."
    "Please choose a supported option: BLAT or BWA."
include: "rules/processing/post_process.rules"

