from os import getenv, getcwd, path
from subprocess import run, PIPE

def import_sample_info(filePath, sampleColName, delim):
    sampleInfo = {}
    with open(filePath, 'r') as info:
        data = info.readlines()
    listData = [row.rstrip().split(delim) for row in data]
    mcols = listData[0]
    samCol = mcols.index(sampleColName)
    samNames = [row.rstrip().split(delim)[samCol] for row in data[1:]]
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

def get_iguide_version(with_hash = False):
    iguide_version_path = getenv("IGUIDE_DIR", "None")
    
    if iguide_version_path in ["None"]:
        raise SystemExit(
          print("\tIGUIDE_DIR cannot be found as an environmental variable.\n"
                "\tCheck to make sure your iGUIDE environment is active,   \n"
                "\tyou may need to restart your environment, update, or    \n"
                "\treinstall iGUIDE with the install.sh script.")
        )
    else:
        iguide_version_path = iguide_version_path + "/.version"
    
    if not path.exists(iguide_version_path):
        raise SystemExit(
          print("\tiGUIDE version cannot be located. Check environmental\n"
                "\tvariables, such as IGUIDE_DIR, otherwise you may want\n"
                "\tto restart your environment, update, or reinstall    \n"
                "\tiGUIDE using the install.sh script.")
        )
    
    iguide_version = open(iguide_version_path, "r").readlines()[0].rstrip()
    
    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
      )
    commit_str = commit_hash.stdout.decode('utf-8').rstrip()
    
    if with_hash:
        return iguide_version + "+" + commit_str
    else:
        return iguide_version

__version__ = get_iguide_version( with_hash = True )
