import sys
import argparse
import subprocess
import iguidelib
from iguidelib.scripts.run import main as Run
from iguidelib.scripts.setup import main as Setup
#from iguidelib.scripts._config import main as Config
#from iguidelib.scripts.list_samples import main as ListSamples

# Goal is to be able to make command line calls in the form:
# iguide <-h/-v/subcommand> <path/to/config.file> <options> -- <snakemake.options>

def main():

    usage_str = "%(prog)s [-h/--help,-v/--version] <subcommand> <path/to/config.file> <options> -- <snakemake.options>"
    description_str = (
        "subcommands:\n"
        "  setup         \tCreate a new config file for a project using local data.\n"
        "  run          \tExecute the iGUIDE pipeline.\n"
        "  config       \tModify or update iGUIDE config files.\n"
        "  list_samples \tMake a list of samples from a directory.\n"
    ).format(version=iguidelib.__version__)

    parser = argparse.ArgumentParser(
        prog="iguide",
        usage=usage_str,
        description=description_str,
        epilog="For more help, see the docs at http://iguide.readthedocs.io.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False
    )

    parser.add_argument("command", help=argparse.SUPPRESS, nargs="?")
    parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s v{}".format(iguidelib.__version__))

    args, remaining = parser.parse_known_args()

    if args.command == "run":
        Run(remaining)
    elif args.command == "setup":
        Setup(remaining)
    elif args.command == "config":
        Config(remaining)
    elif args.command == "list_samples":
        ListSamples(remaining)
    else:
        parser.print_help()
        sys.stderr.write("Unrecognized command.\n")
