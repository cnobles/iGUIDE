# iGUIDE : Improved Genome-wide Unbiased Identification of Double-strand DNA brEaks
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
if ".csv" in config["Sample_Info"]:
    delim = ","
elif ".tsv" in config["Sample_Info"]:
    delim = "\t"
else:
    raise SystemExit("Sample Info file contain extention '.csv' or '.tsv'.")

# Sample information
sampleInfo = import_sample_info(
    config["Sample_Info"], config["Sample_Name_Column"], delim)

SAMPLES=sampleInfo[config["Sample_Name_Column"]]
TYPES=config["Read_Types"]
READS=config["Genomic_Reads"]

# Working paths
RUN = config["Run_Name"]
ROOT_DIR = config["Install_Directory"]
RUN_DIR = config["Install_Directory"] + "/analysis/" + RUN

# Check for directory paths !!! Not sure if this is going to work just yet.
if not os.path.isdir(ROOT_DIR):
    raise SystemExit("Path to iGUIDE is not found. Check configuration file.")

# Target Rules
rule all:
    input: RUN_DIR + "/output/unique_sites." + RUN + ".csv"

# Architecture Rules
include: "rules/workflow_misc/arch.rules"

# Processing Rules
include: "rules/demultiplex/demulti.rules"
include: "rules/sequence_trim/trim.rules"
if (config["UMItags"]):
    include: "rules/sequence_trim/umitag.rules"
    UMIseqs = sampleInfo["barcode2"]
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

