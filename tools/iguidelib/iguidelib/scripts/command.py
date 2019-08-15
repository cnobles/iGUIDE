import sys
import argparse
import subprocess

from iguidelib import __version__
from iguidelib.scripts.run import main as Run
from iguidelib.scripts.setup import main as Setup
#from iguidelib.scripts.config import main as Config
from iguidelib.scripts.list_samples import main as ListSamples
from iguidelib.scripts.evaluate import main as Eval
from iguidelib.scripts.report import main as Report
from iguidelib.scripts.summary import main as Summary
from iguidelib.scripts.clean import main as Clean
from iguidelib.scripts.hints import main as Hints

def main():

    usage_str = "\n  %(prog)s [-h/--help,-v/--version] <subcommand> <path/to/config.file> <options> -- <snakemake.options>"
    description_str = (
        "list of subcommands:\n\n"
        "  primary:\n"
        "    setup        \tCreate a new config file for a project using local data.\n"
        "    run          \tExecute the iGUIDE pipeline.\n\n"
        "  auxiliary:\n"
        "    list_samples \tOutput a list of samples from a project.\n"
        "    eval         \tEvaluate a set or sets of assimilated iGUIDE outputs.\n"
        "    report       \tGenerate a custom report from iGUIDE output files.\n"
        "    summary      \tGenerate a consise summary from iGUIDE output files.\n\n"
        "  additional:\n"
        "    config       \t[inDev] Modify or update iGUIDE config files.\n"
        "    clean        \tCleanup project directory to reduce size. Keeps terminal files.\n"
        "    hints        \tPrint helpful snakemake options to use during processing.\n\n"
    ).format(version=__version__)

    parser = argparse.ArgumentParser(
        prog = "iguide",
        usage = usage_str,
        description = description_str,
        epilog = "For more help, see the docs at http://iguide.readthedocs.io.",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        add_help = False
    )

    parser.add_argument(
        "command", default = "None", help = argparse.SUPPRESS, nargs = "?"
    )
    
    parser.add_argument(
        "-v", "--version", action = "version",
        version = "%(prog)s {}".format(__version__)
    )

    args, remaining = parser.parse_known_args()
    
    sub_cmds = [
        "setup", "run", "eval", "report", 
        "summary", "config", "list_samples", "clean", "hints"
    ]
    
    if not args.command in sub_cmds:
        parser.print_help()
        if not args.command in ['None']:
            sys.stderr.write("  Unrecognized subcommand, '{}'.\n".format(
                args.command
            ))
        sys.exit(1)

    if args.command == "setup":
        Setup(remaining)
    elif args.command == "run":
        Run(remaining)
    elif args.command == "eval":
        Eval(remaining)
    elif args.command == "report":
        Report(remaining)
    elif args.command == "summary":
        Summary(remaining)
    elif args.command == "config":
        raise SystemExit(
            print(
                "\n  'iguide config' subcommand is currently under development."
                "\n  Checkout https://github.com/cnobles/iGUIDE/ for updates"
                "\n  and announcements. Thanks for using iGUIDE!\n"
            )
        )
        #Config(remaining)
    elif args.command == "list_samples":
        ListSamples(remaining)
    elif args.command == "clean":
        Clean(remaining)
    elif args.command == "hints":
        Hints()
    else:
        parser.print_help()
