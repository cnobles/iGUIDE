from os import getenv, getcwd
from subprocess import run, PIPE
from setuptools import setup, find_packages

def get_iguide_version(with_hash = False):
    iguide_version_path = getenv("IGUIDE_DIR", getcwd()) + "/.version"
    iguide_version = open(iguide_version_path, "r").readlines()[0].rstrip()
    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
      )
    commit_str = commit_hash.stdout.decode('utf-8').rstrip()
    if with_hash:
        return iguide_version + "+" + commit_str
    else:
        return iguide_version

setup(
    name = "iguide",
    version = get_iguide_version(),
    packages = find_packages(),
    entry_points = { 'console_scripts': [
        'iguide = iguidelib.scripts.command:main',
        "ig = iguidelib.scripts.command:main"
    ] }
)
