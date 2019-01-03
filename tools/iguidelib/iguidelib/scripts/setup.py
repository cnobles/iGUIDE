import os
import sys
import argparse
from ruamel.yaml import YAML
import subprocess
from pathlib import Path

#from .list_samples import build_sample_list, MissingMatePairError, SampleFormatError
#from sunbeamlib import config
    
def main(argv=sys.argv):
    """Create a new project directory with necessary subdirectories."""

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "Could not determine Conda prefix. Activate your iGUIDE "
            "environment and try this command again.")

    description_str = (
        "Initialize a new iGUIDE project in a given directory, creating "
        "a new config file and (optionally) a sample list.")
    
    parser = argparse.ArgumentParser(
        "setup", description=description_str)
    # TODO
    # Is the -f option needed?
    parser.add_argument(
        "-f", "--force", help="overwrite files if they already exist",
        action="store_true")
    parser.add_argument(
        "--config", help=(
            "name of config file (%(default)s)"),
        default=os.getenv("IGUIDE_DIR", os.getcwd())+"/configs/simulation.config.yml", metavar="FILE")
    parser.add_argument(
        "-i", "--iguide_dir", default=os.getenv("IGUIDE_DIR", os.getcwd()),
        help="Path to iGUIDE installation")
    parser.add_argument(
        "--skip_demultiplexing",action="count",help="Use this if your data is already demultiplexed."
        " (Make sure Demult_Dir is set in config file.)")

    # The remaining args (after --) are passed to Snakemake
    args, remaining = parser.parse_known_args(argv)

    snakefile = Path(args.iguide_dir)/"Snakefile"
    if not snakefile.exists():
        sys.stderr.write(
            "Error: could not find a Snakefile in directory '{}'\n".format(
                args.iguide_dir))
        sys.exit(1)

    yaml=YAML(typ='safe')   # default, if not specfied, is 'rt' (round-trip)
    config=yaml.load(open(args.config, "r"))
    analysis_directory=check_existing(Path("analysis/" + config['Run_Name']))
    
    snakemake_args = ['snakemake', str(analysis_directory),
                      '--snakefile', str(snakefile),
                      '--configfile', str(args.config),
                      '--dir', str(args.iguide_dir)] + remaining
    print("Running: "+" ".join(snakemake_args))

    cmd = subprocess.run(snakemake_args)
    
    if args.skip_demultiplexing:
        # TODO
        # Add use of TYPES argument from config file (i.e. only get I1/I2 if specified)
        try:
            sampleInfo=open(config['Sample_Info'])
        except FileNotFoundError:
            sampleInfo=open(os.getenv("IGUIDE_DIR", os.getcwd())+"/"+config['Sample_Info'])
        sampleList=get_sample_list(sampleInfo)
        demultDir=check_existing(Path(config['Demult_Dir']))
        for sample in sampleList:
            ln_args = ['ln', '-s', str(demultDir)+'/'+sample+'.R1.fastq.gz',
                        str(analysis_directory)+'/processData/'+sample+'.R1.fastq.gz']
            subprocess.run(ln_args)
            ln_args = ['ln', '-s', str(demultDir)+'/'+sample+'.R2.fastq.gz',
                        str(analysis_directory)+'/processData/'+sample+'.R2.fastq.gz']
            subprocess.run(ln_args)
            ln_args = ['ln', '-s', str(demultDir)+'/'+sample+'.I1.fastq.gz',
                        str(analysis_directory)+'/processData/'+sample+'.I1.fastq.gz']
            subprocess.run(ln_args)
            ln_args = ['ln', '-s', str(demultDir)+'/'+sample+'.I2.fastq.gz',
                        str(analysis_directory)+'/processData/'+sample+'.I2.fastq.gz']
            subprocess.run(ln_args)
    else:
        check_existing_fastq(Path(config['R1']))
        check_existing_fastq(Path(config['R2']))
        check_existing_fastq(Path(config['I1']))
        check_existing_fastq(Path(config['I2']))

    
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
