from os import getenv, getcwd, path, chdir
from subprocess import run, PIPE
from pathlib import Path


def import_sample_info(filePath, sampleColName, delim):
    sampleInfo = {}
    with open(filePath, 'r') as info:
        data = info.readlines()
    listData = [row.replace('"', '').rstrip().split(delim) for row in data]
    mcols = listData[0]
    samCol = mcols.index(sampleColName)
    samNames = [row.replace('"', '').rstrip().split(delim)[samCol] for row in data[1:]]
    for m in mcols:
        ind = mcols.index(m)
        vals = []
        for row in listData[1:]:
            vals.append(row[ind])
            colData = dict(zip(samNames, vals))
            sampleInfo[m] = colData
    return sampleInfo


def choose_sequence_data(config_input, sampleInfo):
    if "sampleInfo" in config_input:
        colnam = config_input.split(":")[1]
        if not colnam in sampleInfo:
            raise SystemExit(print("Cannot find ", colnam, "in sampleInfo."))
        seq = sampleInfo[colnam]
    else:
        initial_col = list(sampleInfo)[0]
        samples = list(sampleInfo[initial_col])
        seq = dict(zip(samples, [config_input] * len(samples)))
    return seq


def get_file_path(param, config, root):
    if not str(param) in config:
        raise SystemExit(
            "\n  Cannot locate config parameter: {} \n".format(str(param))
        )
    file_path = Path(config[str(param)])
    if not file_path.exists():
        abs_file_path = Path(root) / file_path
        if not abs_file_path.exists():
            raise SystemExit(
                "\n  Cannot locate file specified by: {} \n".format(config[str(param)])
            )
        else:
            file_path = abs_file_path
    return file_path.absolute()


def get_iguide_version(with_hash = False):
    iguide_path = getenv("IGUIDE_DIR", None)

    if iguide_path is None:
        raise SystemExit(
            print(
                "\n  IGUIDE_DIR cannot be found as an environmental variable."
                "\n  Check to make sure your iGUIDE environment is active,"
                "\n  you may need to restart your environment, update, or"
                "\n  reinstall iGUIDE with the install.sh script.\n"
            )
        )
    else:
        iguide_version_path = iguide_path + "/.version"

    if not path.exists(iguide_version_path):
        raise SystemExit(
            print(
                "\n  iGUIDE version cannot be located. Check environmental"
                "\n  variables, such as IGUIDE_DIR, otherwise you may want"
                "\n  to restart your environment, update, or reinstall"
                "\n  iGUIDE using the install.sh script.\n"
            )
        )

    iguide_version = open(iguide_version_path, "r").readlines()[0].rstrip()

    wd = getcwd()
    chdir(str(iguide_path))

    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
    )

    chdir(wd)

    commit_str = commit_hash.stdout.decode('utf-8').rstrip()

    if with_hash:
        return iguide_version + "+" + commit_str
    else:
        return iguide_version

__version__ = get_iguide_version( with_hash = True )
