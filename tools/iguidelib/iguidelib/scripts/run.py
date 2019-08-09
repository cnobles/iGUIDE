import os
import sys
import argparse
import subprocess

from ruamel.yaml import YAML
from pathlib import Path

def main( argv = sys.argv ):
    """Initiate an iGUIDE project run in Snakemake."""

    try:
        conda_prefix = os.environ.get("CONDA_PREFIX")
    except (KeyError, IndexError):
        raise SystemExit(
            "\n  Could not determine Conda prefix. Activate your iGUIDE "
            "\n  environment and try this command again.\n"
        )

    usage_str = "\n  iguide %(prog)s <path/to/config.file> <options> -- <snakemake.options>"

    description_str = (
        "Initiate the processing of an iGUIDE project givin a configuration "
        "file. Arguments after '--' are passed to Snakemake asis."
    )
    
    parser = argparse.ArgumentParser(
        prog = "run", 
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
        help = "Path to iGUIDE installation."
    )

    # The remaining args (after --) are passed to Snakemake
    args, remaining = parser.parse_known_args(argv)

    snakefile = Path(args.iguide_dir)/"Snakefile"
    if not snakefile.exists():
        sys.stderr.write(
            "\n  Error: could not find a Snakefile in directory '{}'\n".format(
                args.iguide_dir))
        sys.exit(1)

    snakemake_args = ['snakemake',
                      '--snakefile', str(snakefile),
                      '--configfile', str(args.config),
                      '--dir', str(args.iguide_dir)] + remaining
    #print("Running: "+" ".join(snakemake_args))

    cmd = subprocess.run(snakemake_args)
    
    sys.exit(cmd.returncode)
    
        
def check_existing(path, force = False):
    if path.is_dir():
        raise SystemExit(
            "\n  Error: specified file '{}' exists and is a directory.\n".format(path))
    if path.is_file() and not force:
        raise SystemExit(
            "\n  Error: specified file '{}' exists. Use --force to "
            "\n  overwrite.\n".format(path))
    return path
