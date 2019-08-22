def main():

    print(
        "SnakeMake Command Hints:\n\n"
        "  --cores [X]         \tMulticored processing, specified cores to use by X.\n"
        "  --nolock            \tDon't lock iGUIDE directory, for running multiple processes 'iguide run'.\n"
        "  --notemp            \tKeep all temporary files which are defaultly removed during processing.\n"
        "  --keep-going        \tKeep processing even if one job has an error.\n"
        "  -k                  \tShort option for --keep-going.\n"
        "  --latency-wait [X]  \tWait X seconds after a job completes to for output verification, can help with slow servers.\n"
        "  -w [X]              \tShort option for '--latency-wait'.\n"
        "  --restart-times X   \tX is the number of times to restart a job if it fails. Increases 'attempt' each time.\n"
        "  --resources mem_mb=X\tControls the resource limit for 'mem_mb' to help manage pipeline processing.\n"
        "  --rerun-incomplete  \tRe-run all jobs that were not complete before termination.\n"
        "  --ri                \tShort option for '--rerun-incomplete'.\n"
        "  --cluster-config X  \tA JSON or YAML file that defines wildcards used for HPC.\n\n"
        "Remember to pass these options after the '--' flag!\n"
        "Usage:\n"
        "  iguide run <path/to/run.config.yml> -- --nolock --cores 12 --keep-going\n\n"
        "For more help, see the docs at http://iguide.readthedocs.io."
    )
