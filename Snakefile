# iGUIDE : Improved Genome-wide Unbiased Identification of Double-strand DNA brEaks
#
# Author : Christopher Nobles, Ph.D.


import os
import sys
import re
import yaml
import configparser
from iguidelib import import_sample_info, choose_sequence_data

if not config:
    raise SystemExit(
        "No config file specified. Feel free to use the test config as a"
        "template to generate a config file, and specify with --configfile")

# Import sampleInfo
if ".csv" in config["Sample_Info"]:
    delim = ","
elif ".tsv" in config["Sample_Info"]:
    delim = "\t"
else:
    raise SystemExit("Sample Info file needs to contain extention '.csv' or '.tsv'.")

# Sample information
sampleInfo = import_sample_info(
    config["Sample_Info"], config["Sample_Name_Column"], delim)

SAMPLES=sampleInfo[config["Sample_Name_Column"]]
TYPES=config["Read_Types"]
READS=config["Genomic_Reads"]
R1_LEAD=choose_sequence_data(config["R1_Leading_Trim"], sampleInfo)
R1_OVER=choose_sequence_data(config["R1_Overreading_Trim"], sampleInfo)
R2_LEAD=choose_sequence_data(config["R2_Leading_Trim"], sampleInfo)
R2_LEAD_ODN=choose_sequence_data(config["R2_Leading_Trim_ODN"], sampleInfo)
R2_OVER=choose_sequence_data(config["R2_Overreading_Trim"], sampleInfo)

# Working paths
RUN = config["Run_Name"]
ROOT_DIR = ""
try:
    ROOT_DIR = os.environ["IGUIDE_DIR"]
except KeyError:
    raise SystemExit(
        "IGUIDE_DIR environment variable not defined. Are you sure you "
        "activated the iguide conda environment?")
RUN_DIR = ROOT_DIR + "/analysis/" + RUN

# Check for directory paths
if not os.path.isdir(ROOT_DIR):
    raise SystemExit("Path to iGUIDE is not found. Check environmental variables.")

# Target Rules
rule all:
    input: 
      incorp_sites=RUN_DIR + "/output/unique_sites." + RUN + ".csv.gz",
      edit_sites=RUN_DIR + "/output/edited_sites." + RUN + ".rds",
      report=RUN_DIR + "/reports/report." + RUN + ".html",
      stats=RUN_DIR + "/output/stats." + RUN + ".csv"

# Architecture Rules
include: "rules/arch.rules"

# Processing Rules
include: "rules/demulti.rules"
include: "rules/trim.rules"
if (config["UMItags"]):
    include: "rules/umitag.rules"
    UMIseqs = sampleInfo["barcode2"]
include: "rules/filt.rules"
include: "rules/consol.rules"
if (config["Aligner"] == "BLAT" or config["Aligner"] == "blat"):
    include: "rules/align.blat.rules"
    include: "rules/quality.blat.rules"
elif (config["Aligner"] == "BWA" or config["Aligner"] == "bwa"):
    raise SystemExit("BWA aligner not supported yet.")
else:
    "Aligner: " + config["Aligner"] + " not supported."
    "Please choose a supported option: BLAT or BWA."
include: "rules/process.rules"

