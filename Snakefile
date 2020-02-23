# iGUIDE : Improved Genome-wide Unbiased Identification of Double-strand DNA brEaks
#
# Author : Christopher Nobles, Ph.D.


import os
import sys
import re
import yaml
import configparser
from pathlib import Path
from iguidelib import import_sample_info, choose_sequence_data, get_file_path

if not config:
    raise SystemExit(
        "No config file specified. Feel free to use the test config as a"
        "template to generate a config file, and specify with --configfile")

# Working paths
RUN = config["Run_Name"]
ROOT_DIR = ""
try:
    ROOT_DIR = os.environ["IGUIDE_DIR"]
except KeyError:
    raise SystemExit(
        "\n  IGUIDE_DIR environment variable not defined. Are you sure you "
        "\n  activated the iguide conda environment?\n ")
RUN_DIR = ROOT_DIR + "/analysis/" + RUN

# Check for directory paths
if not os.path.isdir(ROOT_DIR):
    raise SystemExit("\n  Path to iGUIDE is not found. Check environmental variables.\n")

# Check for sequence file paths
if not os.path.isdir(config["Seq_Path"]):
    raise SystemExit("\n  Path to sequencing files is not found (Seq_Path). Check your config file.\n")

# Check for config symlink to check proper run directory setup
if not os.path.isfile(RUN_DIR + "/config.yml"):
    raise SystemExit("\n  Path to symbolic config is not present. Check to make sure you've run 'iguide setup' first.\n")

# Check for sampleInfo path
if not "Sample_Info" in config:
    raise SystemExit("\n  Sample_Info parameter missing in config file. Please specify before continuing.\n")
else:
    SAMPLEINFO_PATH = get_file_path("Sample_Info", config, ROOT_DIR)

# Check for suppInfo path
if config["suppFile"]:
    if not "Supplemental_Info" in config:
        raise SystemExit(
            "\n  Supplemental_Info parameter missing in config file."
            "\n  If not including a file, please specify with '.' .\n"
        )
    else:
        if config["Supplemental_Info"] == ".":
            SUPPINFO_PATH = "."
        else:
            SUPPINFO_PATH = get_file_path("Supplemental_Info", config, ROOT_DIR)

# Import sampleInfo
if ".csv" in config["Sample_Info"]:
    delim = ","
elif ".tsv" in config["Sample_Info"]:
    delim = "\t"
else:
    raise SystemExit("\n  Sample Info file needs to contain extention '.csv' or '.tsv'.\n")

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

if not "Alternate_UMI_Method" in config:
    config["Alternate_UMI_Method"] = False
    

# Sample information
sampleInfo = import_sample_info(
    config["Sample_Info"], config["Sample_Name_Column"], delim)

SAMPLES=sampleInfo[config["Sample_Name_Column"]]
READ_TYPES=config["Read_Types"]
READS=config["Genomic_Reads"]
REQ_TYPES=READS[:]

if config["UMItags"] and not config["Alternate_UMI_Method"]: 
    REQ_TYPES.append("I2")

R1_LEAD=choose_sequence_data(config["R1_Leading_Trim"], sampleInfo)
R1_OVER=choose_sequence_data(config["R1_Overreading_Trim"], sampleInfo)
R2_LEAD=choose_sequence_data(config["R2_Leading_Trim"], sampleInfo)
R2_OVER=choose_sequence_data(config["R2_Overreading_Trim"], sampleInfo)

if config["Alternate_UMI_Method"]:
    R1_LEAD_ODN=choose_sequence_data(config["R1_Leading_Trim_ODN"], sampleInfo)
else:
    R2_LEAD_ODN=choose_sequence_data(config["R2_Leading_Trim_ODN"], sampleInfo)


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

if not "qualCtrlMB" in config:
    config["qualCtrlMB"] = 8000
    
if not "assimilateMB" in config:
    config["assimilateMB"] = 4000

if not "evaluateMB" in config:
    config["evaluateMB"] = 4000
    
if not "reportMB" in config:
    config["reportMB"] = 4000

if not "bins" in config:
    config["bins"] = 5

if not "level" in config:
    config["level"] = 300000

if not "readNamePattern" in config:
    config["readNamePattern"] = str("'[\\w\\:\\-\\+]+'")


# Define BINS
BINS = []

for i in range(1, config["bins"] + 1, 1):
    BINS.append("bin" + str(i).zfill(len(str(config["bins"]))))


# Regex constraints on wildcards
wildcard_constraints:
    sample="[\w\-\_]+",
    read="R[12]",
    read_type="[RI][12]",
    req_type="[RI][12]",
    bin="bin[\d]+"

# Target Rules
rule all:
    input: 
      incorp_sites=RUN_DIR + "/output/incorp_sites." + RUN + ".rds",
      report=RUN_DIR + "/reports/report." + RUN + ".html",
      summary=RUN_DIR + "/reports/summary." + RUN + ".txt",
      stats=RUN_DIR + "/reports/runstats." + RUN + ".html"

# Architecture Rules
if (config["Alternate_UMI_Method"]):
    include: "rules/arch.umi_alt_method.rules"
else:
    include: "rules/arch.rules"

# Processing Rules
if (config["skipDemultiplexing"]):
    include: "rules/skip_demulti.rules"
else:
    include: "rules/demulti.rules"
    
include: "rules/binning.rules"

if (config["Alternate_UMI_Method"]):
    include: "rules/trim.umi_alt_method.rules"
else:
    include: "rules/trim.rules"

if (config["UMItags"]):
    if (config["Alternate_UMI_Method"]):
        include: "rules/umitag.umi_alt_method.rules"
    else:
        include: "rules/umitag.rules"
        UMIseqs = sampleInfo["barcode2"]
else:
    include: "rules/umitag_stub.rules"

include: "rules/filt.rules"

if (config["Aligner"] == "BLAT" or config["Aligner"] == "blat"):
    include: "rules/consol.rules"
    include: "rules/align.blat.rules"
    if (config["Alternate_UMI_Method"]):
        include: "rules/quality.blat.umi_alt_method.rules"
    else:
        include: "rules/quality.blat.rules"
elif (config["Aligner"] == "BWA" or config["Aligner"] == "bwa"):
    include: "rules/consol_stub.rules"
    if (config["Alternate_UMI_Method"]):
        include: "rules/align.bwa.umi_alt_method.rules"
    else:
        include: "rules/align.bwa.rules"
    include: "rules/quality.sam.rules"
else:
    raise SystemExit( 
        "\n  Aligner: " + config["Aligner"] + " not currently supported."
        "\n  If you are interested in using the aligner, please contact maintainers."
        "\n  Please choose a supported option: BLAT or BWA.\n"
    )

include: "rules/process.rules"

