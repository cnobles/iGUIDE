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
        default="configs/simulation.config.yml", metavar="FILE")
    parser.add_argument(
        "-i", "--iguide_dir", default=os.getenv("IGUIDE_DIR", os.getcwd()),
        help="Path to iGUIDE installation")

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
    analysis_directory="analysis/" + config['Run_Name']
    
    snakemake_args = ['snakemake', str(analysis_directory), '--snakefile', str(snakefile), '--configfile', str(args.config)] + remaining
    print("Running: "+" ".join(snakemake_args))

    cmd = subprocess.run(snakemake_args)
    
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
