import os
import sys
import argparse
import subprocess

from ruamel.yaml import YAML
from pathlib import Path

def main( argv = sys.argv ):
    """Create a new project directory with necessary subdirectories."""

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your iGUIDE "
            "environment and try this command again.")

    usage_str = "iguide %(prog)s <path/to/config.file> <options> -- <snakemake.options>"
    
    description_str = (
        "Setup a new iGUIDE project given a project configuration file. "
        "Arguments after '--' are passed to Snakemake asis."
    )
    
    parser = argparse.ArgumentParser(
        prog = "setup", 
        usage = usage_str,
        description = description_str
    )

    parser.add_argument(
        "config", 
        help = ("name of config file (%(default)s)"),
        metavar = "CONFIG_FILE"
    )

    parser.add_argument(
        "-i", "--iguide_dir", 
        default = os.getenv("IGUIDE_DIR", os.getcwd()),
        help = "Path to iGUIDE installation")

    parser.add_argument(
        "--skip_demultiplexing", 
        action = "count", 
        help = "Use this option if your data is already demultiplexed."
               " (Make sure Demulti_Dir is set in config file.)")

    # The remaining args (after --) are passed to Snakemake
    args, remaining = parser.parse_known_args(argv)

    snakefile = Path(args.iguide_dir)/"Snakefile"
    
    if not snakefile.exists():
        sys.stderr.write(
            "Error: could not find a Snakefile in directory '{}'\n".format(
                args.iguide_dir))
        sys.exit(1)

    yaml = YAML(typ = 'safe')   # default, if not specfied, is 'rt' (round-trip)
    config = yaml.load(open(args.config, "r"))
    analysis_directory = check_existing(Path("analysis/" + config['Run_Name']))
    read_types = config["Read_Types"]
    
    snakemake_args = ['snakemake', str(analysis_directory),
                      '--snakefile', str(snakefile),
                      '--configfile', str(args.config),
                      '--dir', str(args.iguide_dir)] + remaining
    print("Running: " + " ".join(snakemake_args))

    cmd = subprocess.run(snakemake_args)
    
    if args.skip_demultiplexing:
        try:
            sampleInfo = open(config['Sample_Info'])
        except FileNotFoundError:
            sampleInfo = open(os.getenv("IGUIDE_DIR", os.getcwd()) + "/" + config['Sample_Info'])
        sampleList = get_sample_list(sampleInfo)
        demultiDir = check_existing(Path(config['Demulti_Dir']))
        for sample in sampleList:
            for type in read_types:
                ln_args = [
                    'ln', '-s', str(demultiDir) + '/' + sample + '.' + type + '.fastq.gz',
                    str(analysis_directory) + '/processData/' + sample + '.' + type + '.fastq.gz'
                ]
                subprocess.run(ln_args)
    else:
        for type in read_types: 
            check_existing_fastq(Path(config[type]))
    
    sys.exit(cmd.returncode)
    

def check_existing(path, force=False):
    if path.is_dir():
        raise SystemExit(
            "Error: specified file '{}' exists and is a directory".format(path))
    if path.is_file() and not force:
        raise SystemExit(
            "Error: specified file '{}' exists. Use --force to "
            "overwrite.".format(path))
    return path

def check_existing_fastq(path, force=False):
    if path.is_file() and not force:
        print("Sample file '{}' found.".format(path))
    else:
        print("Warning: specified sample file '{}' does not exist. "
                "Make sure it exists before running iguide run.".format(path))

def get_sample_list(sampleInfo):
    sampleList = []
    for line in sampleInfo:
        items=line.split(',')
        if(items[0]=="sampleName"):
            continue
        sampleList.append(items[0])
    return sampleList
