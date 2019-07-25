# iGUIDE : Improved Genome-wide Unbiased Identification of Double-strand DNA brEaks
#
# Author : Christopher Nobles, Ph.D.


import os
import sys
import re
import yaml
import configparser
from pathlib import Path
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

REQ_READS=READS
if config["UMItags"]: 
    REQ_READS.append("I2")

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

# Check for sequence file paths
if not os.path.isdir(config["Seq_Path"]):
    raise SystemExit("Path to sequencing files is not found (Seq_Path). Check your config file.")

# Check for config symlink to check proper run directory setup
if not os.path.isfile(RUN_DIR + "/config.yml"):
    raise SystemExit("Path to symbolic config is not present. Check to make sure you've run 'iguide setup' first.")

# Default params if not included in config
if not "maxNcount" in config:
    config["maxNcount"] = 1

if not "demultiCores" in config: 
    demulti_cores = snakemake.utils.available_cpu_count()
else:
    demulti_cores = min(
        config["demultiCores"], snakemake.utils.available_cpu_count()
    )

if not "skipDemultiplexing" in config:
    config["skipDemultiplexing"] = False

## Memory and default params
if not "demultiMB" in config:
    config["demultiMB"] = 16000
    
if not "trimMB" in config:
    config["trimMB"] = 4000

if not "filtMB" in config:
    config["filtMB"] = 4000
    
if not "consolMB" in config:
    config["consolMB"] = 4000

if not "alignMB" in config:
    config["alignMB"] = 4000

if not "coupleMB" in config:
    config["coupleMB"] = 4000
    
if not "assimilateMB" in config:
    config["assimilateMB"] = 4000

if not "evaluateMB" in config:
    config["evaluateMB"] = 4000
    
if not "reportMB" in config:
    config["reportMB"] = 4000
  
if not "readNamePattern" in config:
    config["readNamePattern"] = str("'[\\w\\:\\-\\+]+'")

# Target Rules
rule all:
    input: 
      incorp_sites=RUN_DIR + "/output/incorp_sites." + RUN + ".rds",
      report=RUN_DIR + "/reports/report." + RUN + ".html",
      summary=RUN_DIR + "/reports/summary." + RUN + ".txt",
      stats=RUN_DIR + "/reports/runstats." + RUN + ".html"

# Architecture Rules
include: "rules/arch.rules"

# Processing Rules
if (config["skipDemultiplexing"]):
    include: "rules/skip_demulti.rules"
else:
    include: "rules/demulti.rules"
    
include: "rules/trim.rules"

if (config["UMItags"]):
    include: "rules/umitag.rules"
    UMIseqs = sampleInfo["barcode2"]
else:
    include: "rules/umitag_stub.rules"

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

