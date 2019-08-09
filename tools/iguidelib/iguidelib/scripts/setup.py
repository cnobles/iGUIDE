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
            "\n  Could not determine Conda prefix. Activate your iGUIDE "
            "\n  environment and try this command again.\n"
        )

    usage_str = "\n  iguide %(prog)s <path/to/config.file> <options>"
    
    description_str = (
        "Setup a new iGUIDE project given a project configuration file."
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

    # The remaining args will not be used
    args, remaining = parser.parse_known_args(argv)

    # iGUIDE directory
    iguide_directory = Path(args.iguide_dir)
    
    if not iguide_directory.exists():
        sys.stderr.write(
            "Error: could not find iGUIDE directory '{}'.\n".format(
                args.iguide_dir))
        sys.exit(1)
    
    # Load config yaml file
    yaml = YAML(typ = 'safe')
    config = yaml.load(open(args.config, "r"))
    
    analysis_directory = iguide_directory / "analysis" / config['Run_Name']
    
    # Check for existing project directory
    if analysis_directory.exists():
        sys.stderr.write(
            "\n  Error: Project directory currently exists: '{}'.\n".format(
                str(analysis_directory)))
        sys.exit(1)
    else:
        os.makedirs(str(analysis_directory))
    
    # Construct directory tree
    sub_directories = [
        "input_data", "logs", "process_data", "output", "reports"
    ]
    
    for sub_dir in sub_directories:
        os.makedirs(str(analysis_directory / sub_dir))

    # If skipping demultiplexing, create symbolic links to specimen sequencing
    # files.
    read_types = config["Read_Types"] 
     
    if config['skipDemultiplexing']:
        req_read_types = read_types[:]
        req_read_types.remove("I1")
        if not config["UMItags"]: req_read_types.remove("I2")
        try:
            sampleInfo = open(config['Sample_Info'])
        except FileNotFoundError:
            sampleInfo = open(
              os.getenv("IGUIDE_DIR", os.getcwd()) + "/" + config['Sample_Info']
            )
        sampleList = get_sample_list(sampleInfo, config)
        demultiDir = Path(config['Seq_Path'])
        for sample in sampleList:
            for type in req_read_types:
                os.symlink( 
                    str(demultiDir / str(sample + '.' + type + '.fastq.gz')),
                    str(analysis_directory / "input_data" / str(
                        sample + '.' + type + '.fastq.gz'
                        )
                    )
                )
    else:
        for type in read_types: 
            check_existing_fastq(Path(config["Seq_Path"]) / config[type])

    # Create symbolic link to config
    config_path = Path(args.config).absolute()
    
    if config_path.exists():
        os.symlink(str(config_path), str(analysis_directory / "config.yml"))
    else:
        sys.stderr.write(
            "\n  Error: could not locate aboslute path to config file: '{}'.\n".format(
                str(config_path)))
        sys.exit(1)
    
    if analysis_directory.exists():
        print("  '{}' setup has completed.".format(config["Run_Name"]))
    else:
        sys.stderr.write(
            "\n  Error: could not setup project: '{}'.\n".format(
                str(config["Run_Name"])))
        sys.exit(1)
        

def check_existing(path, force=False):
    if path.is_dir():
        raise SystemExit(
            "\n  Error: specified file '{}' exists and is a directory.\n".format(path))
    if path.is_file() and not force:
        raise SystemExit(
            "\n  Error: specified file '{}' exists. Use --force to "
            "\n  overwrite.\n".format(path))
    return path

def check_existing_fastq(path, force=False):
    if path.is_file() and not force:
        print("Sample file '{}' found.".format(path))
    else:
        print("  Warning: specified sample file '{}' does not exist. "
              "Make sure it exists before running iguide run.".format(path))

def get_sample_list(sampleInfo, config):
    sampleList = []
    for line in sampleInfo:
        items=line.split(',')
        if( items[0]==config["Sample_Name_Column"] ):
            continue
        sampleList.append(items[0])
    return sampleList
